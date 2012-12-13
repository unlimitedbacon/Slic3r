package Slic3r::ExPolygon;
use strict;
use warnings;

# an ExPolygon is a polygon with holes

use Boost::Geometry::Utils;
use Math::Geometry::Delaunay 0.08;
use Slic3r::Geometry qw(X Y A B point_in_polygon same_line line_length distance_between_points);
use Slic3r::Geometry::Clipper qw(union_ex JT_MITER);
use List::Util qw(first);

# the constructor accepts an array of polygons 
# or a Math::Clipper ExPolygon (hashref)
sub new {
    my $class = shift;
    my $self;
    if (@_ == 1 && ref $_[0] eq 'HASH') {
        $self = [
            Slic3r::Polygon->new($_[0]{outer}),
            map Slic3r::Polygon->new($_), @{$_[0]{holes}},
        ];
    } else {
        $self = [ map Slic3r::Polygon->new($_), @_ ];
    }
    bless $self, $class;
    $self;
}

sub clone {
    my $self = shift;
    return (ref $self)->new(map $_->clone, @$self);
}

sub contour {
    my $self = shift;
    return $self->[0];
}

sub holes {
    my $self = shift;
    return @$self[1..$#$self];
}

sub lines {
    my $self = shift;
    return map $_->lines, @$self;
}

sub clipper_expolygon {
    my $self = shift;
    return {
        outer => $self->contour,
        holes => [ $self->holes ],
    };
}

sub boost_polygon {
    my $self = shift;
    return Boost::Geometry::Utils::polygon(@$self);
}

sub offset {
    my $self = shift;
    return Slic3r::Geometry::Clipper::offset($self, @_);
}

sub offset_ex {
    my $self = shift;
    return Slic3r::Geometry::Clipper::offset_ex($self, @_);
}

sub safety_offset {
    my $self = shift;
    
    # we're offsetting contour and holes separately
    # because Clipper doesn't return polygons in the same order as 
    # we feed them to it
    
    return (ref $self)->new(
        $self->contour->safety_offset,
        @{ Slic3r::Geometry::Clipper::safety_offset([$self->holes]) },
    );
}

sub noncollapsing_offset_ex {
    my $self = shift;
    my ($distance, @params) = @_;
    
    return $self->offset_ex($distance + 1, @params);
}

sub encloses_point {
    my $self = shift;
    my ($point) = @_;
    return $self->contour->encloses_point($point)
        && (!grep($_->encloses_point($point), $self->holes)
            || grep($_->point_on_segment($point), $self->holes));
}

# A version of encloses_point for use when hole borders do not matter.
# Useful because point_on_segment is slow
sub encloses_point_quick {
    my $self = shift;
    my ($point) = @_;
    return $self->contour->encloses_point($point)
        && !grep($_->encloses_point($point), $self->holes);
}

sub encloses_line {
    my $self = shift;
    my ($line, $tolerance) = @_;
    my $clip = $self->clip_line($line);
    if (!defined $tolerance) {
        # optimization
        return @$clip == 1 && same_line($clip->[0], $line);
    } else {
        return @$clip == 1 && abs(line_length($clip->[0]) - $line->length) < $tolerance;
    }
}

sub point_on_segment {
    my $self = shift;
    my ($point) = @_;
    for (@$self) {
        my $line = $_->point_on_segment($point);
        return $line if $line;
    }
    return undef;
}

sub bounding_box {
    my $self = shift;
    return Slic3r::Geometry::bounding_box($self->contour);
}

sub bounding_box_polygon {
    my $self = shift;
    my @bb = $self->bounding_box;
    return Slic3r::Polygon->new([
        [ $bb[0], $bb[1] ],
        [ $bb[2], $bb[1] ],
        [ $bb[2], $bb[3] ],
        [ $bb[0], $bb[3] ],
    ]);
}

sub bounding_box_center {
    my $self = shift;
    return Slic3r::Geometry::bounding_box_center($self->contour);
}

sub clip_line {
    my $self = shift;
    my ($line) = @_;  # line must be a Slic3r::Line object
    
    return Boost::Geometry::Utils::polygon_linestring_intersection(
        $self->boost_polygon,
        $line->boost_linestring,
    );
}

sub simplify {
    my $self = shift;
    $_->simplify(@_) for @$self;
}

sub scale {
    my $self = shift;
    $_->scale(@_) for @$self;
}

sub translate {
    my $self = shift;
    $_->translate(@_) for @$self;
}

sub rotate {
    my $self = shift;
    $_->rotate(@_) for @$self;
}

sub area {
    my $self = shift;
    my $area = $self->contour->area;
    $area -= $_->area for $self->holes;
    return $area;
}

# Medial axis approximation based on Voronoi diagram.
# Returns a list of polylines and/or polygons.
sub medial_axis {
    my $self = shift;
    my ($width) = @_;
    
    my $expolygon = $self->clone;

    foreach my $polygon (@$expolygon) {
        Slic3r::Geometry::polyline_remove_short_segments($polygon, $width / 2);
        # minimally subdivide polygon segments so that the
        # triangulation captures important shape features, like corners 
        $polygon->splitdivide($width/2);
    }
    
    my $tri = Math::Geometry::Delaunay->new();
    $tri->doVoronoi(1);
    $tri->doEdges(1);
    $tri->addPolygon($expolygon->contour);
    $tri->addHole($_) for $expolygon->holes;
    my ($topo, $vtopo) = $tri->triangulate();

    # Move vertices in Voronoi diagram ($vtopo) to make them 
    # work much better as a medial axis approximation.

    Math::Geometry::Delaunay::mic_adjust($topo, $vtopo);

    # Remove references to ray edges from all nodes.

    foreach my $node (@{$vtopo->{nodes}}) {
        # remove refs to ray edges from nodes
        # vector [0,0] means it's not a ray, so keep it
        @{$node->{edges}} = grep !$_->{vector}->[0] && !$_->{vector}->[1], @{$node->{edges}};
    }
    
    # Get the subset of nodes that are not too close to the polygon boundary.

    my @vnodes;
    my @new_nodes;
    my @adjust_queue = (0 .. $#{$vtopo->{nodes}});

    while (scalar @adjust_queue) {
        my $node_index = shift @adjust_queue;
        my $node = $vtopo->{nodes}->[$node_index];

        if (!defined $node->{radius}) {

            # Distance from any node in the Voronoi diagram to the nearest point on the
            # polygon can be approximated. For each edge emanating from the node, look
            # up the corresponding edge in the Delaunay triangulation. The distance from
            # the node point to either end of the Delaunay edge is roughly the radius
            # of the maximum inscribed circle to the polygon at the node.
            # (The approximate radius is greater than the true radius.)

            for (my $i = $#{$node->{edges}}; $i > -1; $i--) {
                $node->{radius} = distance_between_points($node->{point},$topo->{edges}->[$node->{edges}->[$i]->{index}]->{nodes}->[0]->{point});
                if ($node->{radius} < $width / 2) {
                    my $edge = splice(@{$node->{edges}}, $i, 1);
                }
            }
        } else {

            # If we ran mic_adjust(), each Voronoi node has been shifted
            # and given a radius value for the maximum inscribed circle.
            # Trimming based on this radius should be more accurate than 
            # the check above, which should be a fallback in case radius 
            # could not be determined for some reason.

            if ($node->{radius} < $width / 2) {

                for (my $i = $#{$node->{edges}}; $i > -1; $i--) {

                    my ($other_node_index, $this_node_index) = $node == $node->{edges}->[$i]->{nodes}->[0] ? (1, 0) : (0, 1);
                    my $other_node = $node->{edges}->[$i]->{nodes}->[$other_node_index];

                    if ($other_node->{radius}  > $width / 2) {

                        # All nodes directly connected to this one should maybe 
                        # reevaluate their situation after this one has been processed,
                        # so add them to the end of the list again.
                        for (my $j = $#{$node->{edges}}; $j > -1; $j--) {
                            my $other_node_index = $node == $node->{edges}->[$j]->{nodes}->[0] ? 1 : 0;
                            # convert to the index in $vtopo->{nodes} list
                            $other_node_index = $node->{edges}->[$j]->{nodes}->[$other_node_index]->{index};
                            push @adjust_queue, $other_node_index;
                        }

                        # clone
                        my $this_node = {
                                          index  => scalar(@{$vtopo->{nodes}}),
                                          point  => [@{$node->{point}}],
                                          radius => $node->{radius},
                                          edges  => [],
                                          elements => [],
                                          segments => [],
                                          marker => undef,
                                          attributes => [],
                                          };

                        # attach current edge to new node, detach from old
                        my $edge = $node->{edges}->[$i];
                        splice @{$edge->{nodes}}, $this_node_index, 1, $this_node;
                        @{$this_node->{edges}} = splice @{$node->{edges}}, $i, 1;

                        # In order to have the option to retain some nodes with 
                        # too-small radius in a later step, add an edge that 
                        # links the the new node to the old one.
                        push @{$vtopo->{edges}}, {
                                                  nodes    => [$node,$this_node],
                                                  index    => scalar(@{$vtopo->{edges}}),
                                                  vector   => [0,0],
                                                  elements => [],
                                                  };
                        push @{$node->{edges}}, $vtopo->{edges}->[-1];
                        push @{$this_node->{edges}}, $vtopo->{edges}->[-1];
                        
                        push @new_nodes, $this_node;
                        push @{$vtopo->{nodes}}, $this_node;
                        
                        # interpolate to where radius == width
                        my $factor = (($width / 2) - $this_node->{radius}) / ($other_node->{radius} - $this_node->{radius});
                        $this_node->{point} = [$this_node->{point}->[0] + $factor * ($other_node->{point}->[0] - $this_node->{point}->[0]),
                                               $this_node->{point}->[1] + $factor * ($other_node->{point}->[1] - $this_node->{point}->[1])];
                        $this_node->{radius} = $width / 2;
                    }
                }
            }
        }
    }

    # Remove any nodes that still have too small 
    # of a radius after the adjustments above.

    # Though maybe we want expand this to preserve some 
    # in some cases - like if we want to make thin walls 
    # of extrusion width even when model walls are thinner than that.

    foreach my $node (@{$vtopo->{nodes}}) {
        if ($node->{radius} < $width / 2) {
            for (my $i = $#{$node->{edges}}; $i > -1; $i--) {
                my $other = $node == $node->{edges}->[$i]->{nodes}->[0]
                            ? $node->{edges}->[$i]->{nodes}->[1]
                            : $node->{edges}->[$i]->{nodes}->[0];
                @{$other->{edges}} = grep $_ != $node->{edges}->[$i] ,@{$other->{edges}};
            }
            $node->{edges} = [];
        }
    }

    push @vnodes, grep @{$_->{edges}}, @{$vtopo->{nodes}};

    # Branch with short twigs filter
    # Often a path ends in a legitimate fork, with two very
    # short twigs coming out of a branch node. We want to
    # convert              /    to    __---    ,
    #                 \   /           \
    #                  \ /             \
    #                   |               |
    #                   |               |
    # to eliminate the branch, for a continuous toolpath.

    foreach my $branch_node (grep @{$_->{edges}} > 2, @vnodes) {
        my @twig_tips;
        foreach my $edge (@{$branch_node->{edges}}) {
            # Walk out on each edge to find short (but possibly multi-edge) twigs.
            my @prune;
            my $tip = +(grep $_ != $branch_node && @{$_->{edges}} > 0,  @{$edge->{nodes}})[0];
            next if !$tip;
            my $next_edge = $edge;
            my $dist = distance_between_points($tip->{point},$branch_node->{point});
            push @prune, $tip;
            while (@{$tip->{edges}} != 1 && $dist < $width) {
                $next_edge = +(grep $_ != $next_edge, @{$tip->{edges}})[0];
                $tip = +(grep $_ != $tip && @{$_->{edges}} > 0, @{$next_edge->{nodes}})[0];
                if (!$tip) {last;}
                push @prune, $tip;
                $dist = distance_between_points($tip->{point},$branch_node->{point});
            }
            push @twig_tips, [@prune,$dist] if (@{$prune[-1]->{edges}} == 1 && $dist < $width);
        }
        # Branch node will connect to closest tip,
        # then closest will connect to next closest.
        @twig_tips = sort {$a->[-1] <=> $b->[-1]} @twig_tips;
        if (@twig_tips > 1) {
            my $dist0  = pop @{$twig_tips[0]};
            my $dist1  = pop @{$twig_tips[1]};

            # Remove from branch node the edge ref that lead to tip1.
            foreach my $e (@{$twig_tips[1]->[0]->{edges}}) {
                @{$branch_node->{edges}} = grep $_ != $e, @{$branch_node->{edges}};
            }

            # Make tip1's only edge connect to tip0.
            my $tip0 = pop @{$twig_tips[0]};
            my $tip1 = pop @{$twig_tips[1]};
            $tip1->{edges}->[0]->{nodes} = [$tip1,$tip0];
            push @{$tip0->{edges}}, $tip1->{edges}->[0];

            # Prune any nodes that were between tip1 and the branch node.
            $_->{edges} = [] for @{$twig_tips[1]};
        }        
    }

    @vnodes = grep @{$_->{edges}} > 0, @vnodes;

    # All nodes where more than two edges meet are branch nodes.
    my @branch_start_nodes = grep @{$_->{edges}} > 2, @vnodes;

    # If no branch nodes, we're dealing with a line or loop.
    # If a line, nodes at ends will have only one edge reference. 
    if (@branch_start_nodes == 0) {
        @branch_start_nodes = grep @{$_->{edges}} == 1, @vnodes;
    }

    # Otherwise, it's a loop - any node can be the start node.
    if (@branch_start_nodes == 0) {
        push @branch_start_nodes, $vnodes[0];
    }

    my @polyedges = ();
    my @end_edges = ();
    
    # Walk the cross referenced nodes and edges to build up polyline-like node lists.
    foreach my $start_node (@branch_start_nodes) {
        foreach my $start_edge (@{$start_node->{edges}}) {
            # don't go back over path already traveled
            next if first {$_ == $start_edge} @end_edges;
            my $this_node = $start_node;
            push @polyedges, [];
            push @{$polyedges[-1]}, $this_node;
            my $this_edge = $start_edge;
            #step along nodes: next node is the node on current edge that isn't this one
            while ($this_node = +(grep $_ != $this_node, @{$this_edge->{nodes}})[0]) {
                # stop at point too close to polygon
                last if (@{$this_node->{edges}} == 0);
                # otherwise, always add the point - duplicate start and end lets us detect polygons
                push @{$polyedges[-1]}, $this_node;
                # stop at a branch node or dead end, and remember the edge so we don't backtrack
                if (@{$this_node->{edges}} > 2 || @{$this_node->{edges}} == 1 ) {
                    push @end_edges, $this_edge;
                    last;
                }
                # step to next edge
                $this_edge = +(grep $_ != $this_edge, @{$this_node->{edges}})[0];
                # stop if we've looped around to start
                last if $this_edge == $start_edge;
            }
        }
    }

    # Now combine chains of edges into longer chains where their ends meet,
    # deciding which two chains to link at branch node sites.

    my @polylines = ();

    # Sort by length, where array length is rough proxy for edge length sum.
    @polyedges = sort {@{$a} <=> @{$b}} @polyedges;

    # Link polyedges with common end points.
    for (my $i = $#polyedges; $i > 0; $i--) {
        # polygons
        if ($polyedges[$i]->[0] == $polyedges[$i]->[@{$polyedges[$i]} - 1]) {
            push @polylines, splice(@polyedges, $i, 1);
            next;
        }
        # polylines
        else {

            # We take a longer polyline from the end of the list
            # and see if it links up with any of the shorter
            # polylines that come before it in the list.
            # If so, we splice the longer polyline out of the list
            # and add it's points to the shorter polyline.
            # Having that first splice choose the shorter next path
            # is meant to fill out local features while in the neighborhood
            # instead of always linking longest paths, which might
            # require revisting a lot of separate regions of local features
            # that were passed by, requiring more rapid traversals.
            # The paths may need further sorting after this linking to achieve
            # this. The point here is just to link up the path structure that 
            # will enable that.
            # ... sort of works - but there's probably a better approach
            
            my $this  = $polyedges[$i];

            for (my $j = 0; $j < $i ; $j++) {
                my $other = $polyedges[$j];
                # all the cases of ends matching up
                if ($this->[@{$this} - 1] == $other->[0]) {
                    shift @{$other};
                    @{$other} = (@{splice(@polyedges, $i, 1)}, @{$other});
                    last;
                } elsif ($this->[0] == $other->[@{$other} - 1]) {
                    shift @{$this};
                    @{$other} = (@{$other}, @{splice(@polyedges, $i, 1)});
                    last;
                } elsif ($this->[0] == $other->[0]) {
                    shift @{$this};
                    @{$other} = ((reverse @{$other}), @{splice(@polyedges, $i, 1)});
                    last;
                } elsif ($this->[@{$this} - 1] == $other->[@{$other} - 1]) {
                    pop @{$other};
                    @{$other} = (@{splice(@polyedges, $i ,1)}, (reverse @{$other}));
                    last;
                }
            }
        }
    }

    push @polylines, @polyedges;

    my @result = ();
    foreach my $polyline (@polylines) {
        next unless @$polyline >= 2;
        
        @$polyline = map $_->{point}, @$polyline;
                        
        if (Slic3r::Geometry::same_point($polyline->[0], $polyline->[-1])) {
            next if @$polyline == 2;
            push @result, Slic3r::Polygon->new(@$polyline[0..$#$polyline-1]);
        } else {
            push @result, Slic3r::Polyline->new($polyline);
        }
    }
    
    return @result;
}

package Slic3r::ExPolygon::Collection;
use Moo;
use Slic3r::Geometry qw(X1 Y1);

has 'expolygons' => (is => 'ro', default => sub { [] });

sub clone {
    my $self = shift;
    return (ref $self)->new(
        expolygons => [ map $_->clone, @{$self->expolygons} ],
    );
}

sub align_to_origin {
    my $self = shift;
    
    my @bb = Slic3r::Geometry::bounding_box([ map @$_, map @$_, @{$self->expolygons} ]);
    $_->translate(-$bb[X1], -$bb[Y1]) for @{$self->expolygons};
}

sub rotate {
    my $self = shift;
    $_->rotate(@_) for @{$self->expolygons};
}

sub size {
    my $self = shift;
    return [ Slic3r::Geometry::size_2D([ map @$_, map @$_, @{$self->expolygons} ]) ];
}

1;
