"""

Iso-Delaunay

"""


from flatsurf import *
H = HyperbolicPlane(AA)
from surface_dynamics import *
from flatsurf.geometry.euclidean import ccw
# our first set of import. We may want to change it later


class FlatTriangulation:
    r"""
    EXAMPLES :

        sage: from flatsurf import translation_surfaces
        sage: S = translation_surfaces.arnoux_yoccoz(3)
        sage: FlatTriangulation(S.triangulate())
        Flat Triangulation of Translation Surface in H_3(2^2) built from 12 triangles

    The surface must be triangulated.

        sage: S = MutableOrientedSimilaritySurface(AA)
        sage: p0 = polygons.regular_ngon(4)
        sage: S.add_polygon(p0)
        sage: S.glue((0,0),(0,2))
        sage: S.glue((0,1),(0,3))
        sage: S.set_immutable()
        sage: FlatTriangulation(S)
        ValueError: The surface is not triangulated
        sage: FlatTriangulation(S.triangulate())
        Flat Triangulation of Translation Surface in H_1(0) built from a square
    
    """
    def __init__(self,M):
        if not M.is_triangulated():
            raise ValueError("The surface is not triangulated")
        self._surface = M
        edges = M.edges()
        m = len(edges)
        fp = [-1]*m
        ep = [2*(i//2)+(1- i%2) for i in range(m)]
        self._vectors = [] # the list of all the vectors of the edges : it contains the information of the holonomy and of the triangle (from the initial triangulation of M) it belongs
        labelling = {}
        for (a,b) in edges:
            if (a,b) not in labelling:
                (c,d) = M.opposite_edge(a,b)
                labelling[a,b]=len(labelling)
                self._vectors.append(M.tangent_vector(a,M.polygon(a).vertex(b),M.polygon(a).edge(b)))
                labelling[c,d]=len(labelling)
                self._vectors.append(M.tangent_vector(c,M.polygon(c).vertex(d),M.polygon(c).edge(d)))
        for lab in M.labels():
            for i in range(3):
                fp[labelling[lab, i]]=labelling[lab, (i+1)%3]

        # NOTE: on recent versions of surface_dynamics, ep is forced to be
        # 0 <-> 1, 2 <-> 3, etc and is not anymore an argument of the constructor
        try:
            self._triangulation = FatGraph(ep=ep, fp=fp) # the graph given by the triangulation
        except TypeError:
            self._triangulation = FatGraph(fp=fp) # the graph given by the triangulation

    def _check(self):
        self._triangulation._check()

        fp = self._triangulation.face_permutation(copy=False)
        n = len(self._vectors)
        for e1 in range(0, n, 2*n):
            v1 = self._vectors[e1].vector()
            if v1 != -self._vectors[e1 + 1].vector():
                raise ValueError('edge ({}, {}) have vectors with non-opposite holonomies'.format(e1, e1+1))

            e2 = fp[e1]
            v2 = self._vectors[e2].vector()
            e3 = fp[e2]
            v3 = self._vectors[e3].vector()
            if not (v1 + v2 + v3).is_zero() or ccw(v1, v2) <= 0 or ccw(v2, v3) <= 0 or ccw(v3, v1) <= 0:
                raise ValueError('({}, {}, {}) is not a valid triangle'.format(e1, e2, e3))

    def __repr__(self):
        return 'Flat {}'.format(self._surface)

    def copy(self):
        nt = FlatTriangulation(self._surface)
        nt._vectors = self._vectors[:]
        nt._triangulation = self._triangulation.copy()
        return nt

    def is_isomorphic(self, other, certificate=False):
        r"""
        Test whether ``self`` and ``other`` defines the same triangulation of the underlying surface.

        OUTPUT: if ``certificate=False`` then the output is a boolean. If ``certificate=True`` then
        the output is a pair ``(answer, relabelling)`` where ``answer`` is a boolean and ``relabelling``
        is ``None`` if ``answer=False`` and is a permutation of labels if ``answer=True``.

        EXAMPLES::

            sage: from flatsurf import translation_surfaces
            sage: s = translation_surfaces.arnoux_yoccoz(3)
            sage: t0 = FlatTriangulation(s)
            sage: t = t0.copy()
            sage: t0.is_isomorphic(t)
            True
            sage: t0.is_isomorphic(t, certificate=True)
            (True, [0, 1, 2, 3, ..., 35])

        Now flip an edge twice so that the labels get flipped::

            sage: e = choice([e for e in range(2 * t.triangulation.num_edges()) if t.is_flipable(e)])
            sage: t.flip(e)
            sage: t0.is_isomorphic(t)
            False
            sage: t.flip(e)
            sage: t0.is_isomorphic(t)
            True
            sage: ans, certificate = t0.is_isomorphic(t, certificate=True)
            sage: ans
            True
            sage: e != certificate[e]
            True
        """
        if not isinstance(other, FlatTriangulation):
            raise TypeError('can not compare t of type {} with a FlatTriangulation'.format(type(t).__name__))
        if self._surface is not t._surface:
            raise ValueError('can not compare FlatTriangulation on different underlying surfaces')

        # find the first vector of self in other
        n = len(self._vectors)
        u0 = self._vectors[0]
        i = 0
        while i < n:
            v = other._vectors[i]
            if u0 == v:
                break
            i += 1
        if i == n:
            return (False, None) if certificate else False

        # does simultaneous DFS
        sfp = self._triangulation.face_permutation(copy=False)
        ofp = other._triangulation.face_permutation(copy=False)
        relabelling = [None] * n
        relabelling[0] = i
        relabelling[1] = i ^ 1
        todo = [0, 1]
        while todo:
            i = todo.pop()
            ii = relabelling[i]
            assert ii is not None

            for _ in range(3):
                if relabelling[i] is None:
                    if self._vectors[i] != other._vectors[ii]:
                        return (False, None) if certificate else False
                    relabelling[i] = ii
                    relabelling[i ^ 1] = ii ^ 1
                    todo.append(i ^ 1)
                elif relabelling[i] != ii:
                    return (False, None) if certificate else False
                i = sfp[i]
                ii = ofp[ii]

        return (True, relabelling) if certificate else True

    def graph(self):
        r"""
        Return the FatGraph (from surfaces_dynamics) that describes the triangulations
        """
        return self._triangulation.copy()

    def vectors(self):
        r"""
        Return the list of the saddle connections of the Flat Triangulation
        """
        return self._vectors[:]

    def is_flipable(self,e1):
        r"""
        Test whether the edge can be flipped to form a new triangulation

        EXAMPLES :
            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.arnoux_yoccoz(3)
            sage: t = FlatTriangulation(S.triangulate())
            sage: t.is_flipable(0)
            True
            sage: t.flip(0)
            sage:t.is_flipable(4)
            False
        """
        if e1%2 == 0:
            e2 = e1+1
        else : e2 = e1-1
        fp = self._triangulation.face_permutation(copy = False)
        f1 = fp[e1]
        f2 = fp[f1]
        g1 = fp[e2]
        g2 = fp[g1]
        return (ccw(self._vectors[f2].vector(),self._vectors[g1].vector()) > 0) and (ccw(self._vectors[g2].vector(),self._vectors[f1].vector()) > 0)

    def flip(self, e1):
        r"""
        Return the FlatTriangulation obtained after the flip of the edge ``e1``.

        EXAMPLES :
        
            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.arnoux_yoccoz(3)
            sage: t = FlatTriangulation(S.triangulate())
            sage: t.is_flipable(0)
            sage: t.flip(0)
            sage: t
            Flat Triangulation of Translation Surface in H_3(2^2) built from 12 triangles

        The edge must be flipable.
        
            sage: t.flip(4)
            ValueError: This edge is not flipable
        """
        
        M = self._surface
        if not self.is_flipable(e1):
             raise ValueError("This edge is not flipable")
        # modification of the graph
        if e1%2 == 0:
            e2 = e1+1
        else : e2 = e1-1
        fp = self._triangulation.face_permutation(copy=True)
        f1 = fp[e1]
        f2 = fp[f1]
        g1 = fp[e2]
        g2 = fp[g1]
        fp[e1] = f2
        fp[f2] = g1
        fp[g1] = e1
        fp[e2] = g2
        fp[g2] = f1
        fp[f1] = e2
        # TODO: in order to not mess up the internal structure of FatGraph we build a new
        # FatGraph... which is not optimal
        ep = [2*(i//2)+(1- i%2) for i in range(len(M.edges()))]
        try:
            self._triangulation = FatGraph(ep=ep, fp=fp)
        except TypeError:
            self._triangulation = FatGraph(fp=fp)

        # modification of the vectors
        hol2 = self._vectors[g1].vector()+self._vectors[f2].vector()
        hol1 = -hol2
        
        lab_pol = self._vectors[g2].polygon_label()
        lab_ver = self._vectors[g2].vertex()
        lab_edg = (lab_ver-1)%3 # we will check if the edge of the polygon lab_pol is ccw of the vector or not
        vtest = -M.polygon(lab_pol).edge(lab_edg)
        while ccw(vtest,hol1)>0:
            (nf,ne) = M.opposite_edge(lab_pol,lab_edg)
            lab_pol = nf
            lab_edg = (ne-1)%3
            lab_ver = ne
            vtest = -M.polygon(lab_pol).edge(lab_edg)
        pol1 = lab_pol
        ver1 = lab_ver


        lab_pol = self._vectors[f2].polygon_label()
        lab_ver = self._vectors[f2].vertex()
        lab_edg = (lab_ver-1)%3 # we will check if the edge of the polygon lab_pol is ccw of the vector or not
        vtest = -M.polygon(lab_pol).edge(lab_edg)
        while ccw(vtest,hol2)>0:
            (nf,ne) = M.opposite_edge(lab_pol,lab_edg)
            lab_pol = nf
            lab_edg = (ne-1)%3
            lab_ver = ne
            vtest = -M.polygon(lab_pol).edge(lab_edg)
        pol2 = lab_pol
        ver2 = lab_ver

        self._vectors[e1] = M.tangent_vector(pol1, M.polygon(pol1).vertex(ver1), hol1)
        self._vectors[e2] = M.tangent_vector(pol2, M.polygon(pol2).vertex(ver2), hol2)

    
    def geodesic(self, e1): 
        r"""
        Return the hyperbolic geodesic associate to the flip of the edge with label ``e1`` in the FlatTriangulation.
        This function is used to compute the IsoDelaunayCell class.

        EXAMPLES :
        
            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.arnoux_yoccoz(3)
            sage: t = FlatTriangulation(S.triangulate())
            sage: t.geodesic(0)
            {-0.0797431060641191?*x - 0.08737802538415273? = 0}

        """
        
        if not self.is_flipable(e1):
            return None
        if e1%2 == 0:
            e2 = e1+1
        else : e2 = e1-1
        fp = self._triangulation.face_permutation(copy = False)
        f1 = fp[e1]
        f2 = fp[f1]
        g1 = fp[e2]
        g2 = fp[g1]
        x1 = self._vectors[g1].vector()[0]
        y1 = self._vectors[g1].vector()[1]
        x2 = self._vectors[e1].vector()[0]
        y2 = self._vectors[e1].vector()[1]
        V = -self._vectors[f2].vector()
        x3 = V[0]
        y3 = V[1]
        a = x1*y2*y3*(y3-y2) + x2*y1*y3*(y1-y3) + x3*y1*y2*(y2-y1)
        b = x1*y1*(x2*y3-x3*y2) + x2*y2*(x3*y1-x1*y3) + x3*y3*(x1*y2-x2*y1)
        c = x1*x2*y3*(x1-x2) + x2*x3*y1*(x2-x3) + x1*x3*y2*(x3-x1)
        return H.geodesic(a,2*b,c=c,model = "half_plane")
    
    def plot(self, *args, **kwds):
        res = Graphics()
        for v in self._vectors[::2]:
            traj = v.straight_line_trajectory()
            traj.flow(100)
            res = res + traj.plot(*args, **kwds)
        return res





class IsoDelaunayCell:
    r"""
    A cell of the hyperbolic plane associated to a FlatTriangulation.

    EXAMPLES :

        sage: from flatsurf import translation_surfaces
        sage: S = translation_surfaces.arnoux_yoccoz(3)
        sage: t = FlatTriangulation(S.triangulate())
        sage: IsoDelaunayCell(t)
        Iso-Delaunay Cell of a Flat Triangulation of Translation Surface in H_3(2^2) built from 12 triangles
    
    """
    def __init__(self,t):
        self._triangulation = t.copy()
        Passage1 = True
        n = len(t._vectors)
        fp = t._triangulation.face_permutation(copy = False)
        cor = {}
        for i in range(0,n,2):
            g = t.geodesic(i)
            if g in cor:
                cor[g].append(i)
            else:
                cor[g] = [i]
            if g != None:
                if Passage1:
                    Sol = g.left_half_space()
                    Passage1 = False
                else:
                    Sol = H.intersection(Sol, g.left_half_space())
        self._polygon = Sol
        e = [s.geodesic() for s in Sol.edges()]
        cor_final = {}
        for g in e:
            cor_final[g] = cor[g]
        self._correspondance = cor_final

    def __repr__(self):
        return 'Iso-Delaunay Cell of a {}'.format(self._triangulation)

    def __eq__(self, other):
        if not isinstance(other, IsoDelaunayCell):
            raise ValueError('can not compare object of type {} with IsoDelaunayCell'.format(type(other).__name__))
        if not self._triangulation._surface == other._triangulation._surface:
            raise ValueError('can not compare IsoDelaunayCell from different surface')
        return self._polygon == other._polygon

    def edges(self):
        r"""
        Return the list of geodesic that defines the cell in ccw order.
        This geodesics are oriented !
        """
        e = self._polygon.edges()
        return [elt.geodesic() for elt in e]

    def adjacent(self,other):
        r"""
        Test whether two IsoDelaunayCell are adjacents meaning that they share a common edge.

        EXAMPLES :

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.arnoux_yoccoz(3)
            sage: t = FlatTriangulation(S.triangulate())
            sage: R1 = IsoDelaunayCell(t)
            sage: t.flip(2)
            sage: R2 = IsoDelaunayCell(t)
            sage: R1.is_adjacent(R2)
            True
            sage: t.flip(6)
            sage: t.flip(22)
            sage: R3 = IsoDelaunayCell(t)
            sage: R1.adjacent(R3)
            False
        """
        
        if not isinstance(other, IsoDelaunayCell):
            raise ValueError('can not compare adjacency of object of type {} with IsoDelaunayCell'.format(type(other).__name__))
        if not self._triangulation._surface == other._triangulation._surface:
            raise ValueError('can not compare adjacency of IsoDelaunayCell from different surface')
        if self == other:
            return False
        p1 = self._polygon
        p2 = other._polygon
        i = p1.intersection(p2)
        if i.is_empty():
            return False
        elif i.dimension() == 0:
            return False
        elif i.dimension() == 1:
            return True
        else:
            return ValueError('this case should not occur.')

    def area(self, numerical = True):
        return self._polygon.area(numerical)
    
    def plot(self, *args, **kwds):
        return self._polygon.plot(*args, **kwds)


class IsoDelaunayTessellation:
    r"""
    Tesselation of the hyperbolic plane from a certain surface.
    Will contain the information of every explored IsoDelaunayCell.

    EXAMPLES :

        sage: from flatsurf import translation_surfaces
        sage: S = translation_surfaces.arnoux_yoccoz(3)
        sage: t = FlatTriangulation(S.triangulate())
        sage: R = IsoDelaunayCell(t)
        sage: IsoDelaunayTessellation(R)
        Iso-Delaunay tesselation of the hyperbolic plane from 1 iso-Delaunay cells of a Triangulation of Translation Surface in H_3(2^2) built from 12 triangles
    """
    def __init__(self,c):
        self._explored = {} #dictionnary of the cell already explored
        self._explored[0] = c
        k = len(c._correspondance)
        ep = [2*(i//2)+(1- i%2) for i in range(2*k)]
        fp = [-1]*(2*k)
        geo = c.edges()
        self._edge_to_vertex = [None]*(2*k)
        self._edge_to_geodesic = [None]*(2*k)
        vp = [-1]*(2*k)
        for i in range(k):
            geo[i] = geo[i].geodesic()
            fp[2*i] = (2*i+3)%(2*k)
            fp[(2*i+3)%(2*k)] = 2*i
            self._edge_to_vertex[2*i] = 0
            self._edge_to_vertex[(2*i+3)%(2*k)] = 1
            self._edge_to_geodesic[2*i] = geo[i].geodesic()
            vp[2*i] = (2*i+2)%(2*k)
            vp[2*i+1] = (2*i+3)%(2*k)
        try:
            self._graph = FatGraph(ep=ep, fp=fp, mutable = True) # the dual-graph of the cell already explored
        except TypeError:
            self._graph = FatGraph(fp=fp, mutable = True)
        self._boundary = [1] # list of the boundary vertex (in the FatGraph, it correponds to a vertex but it's a Delaunay cell)

    def __repr__(self):
        return 'Iso-Delaunay tesselation of the hyperbolic plane from {} iso-Delaunay cells of a {}'.format(len(self._explored),self._explored[0]._triangulation._surface)

    def _check(self):
        self._graph._check()
        if self._graph.euler_characteristic() != 2:
            raise ValueError('The graph should be planar')
        n = self._graph.num_vertices()
        if n > len(self._explored)+len(self._boundary):
            raise ValueError('There is vertex that are neither explored and unexplored')
        elif n < len(self._explored)+len(self._boundary):
            raise ValueError('There is a vertex that is simultanely explored and unexplored')
        for elt in self._boundary:
            if elt in self._explored.keys():
                raise ValueError('There is a vertex that is simultanely explored and unexplored')
        m = len(self._graph.face_permutation())
        if m != len(self._edge_to_vertex) or m != len(self._edge_to_geodesic):
            raise ValueError('There is an edge that is not assigned')
        for e in self.explorable_edges()[0]:
            v = self._edge_to_vertex[e]
            geo = self._edge_to_geodesic[e]
            oc = self._explored[v]
            try : oc._correspondance[geo]
            except KeyError : raise ValueError('There is an edge whose correspondance is not assigned')

    def explore(self, e):
        r"""
        Add the cell obtaining by exploring accross the semi-edge ``e``

        EXAMPLES :

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.arnoux_yoccoz(3)
            sage: t = FlatTriangulation(S.triangulate())
            sage: R = IsoDelaunayCell(t)
            sage: T = IsoDelaunayTessellation(R)
            sage: T.explore(0)
            sage: T
            Iso-Delaunay tesselation of the hyperbolic plane from 2 iso-Delaunay cells of a Triangulation of Translation Surface in H_3(2^2) built from 12 triangles

        The edge must be explorable.

            sage: T.explore(1)
            ValueError: The vertex linked to this edge has already been explored
            sage: T.explore(5)
            ValueError: This edge come from an unexplored vertex
        """
        if self._edge_to_geodesic[e] == None:
            raise ValueError('This edge come from an unexplored vertex')
        if e%2 == 0:
            e2 = e+1
        else : e2 = e-1
        if self._edge_to_geodesic[e2] != None:
            raise ValueError("The vertex linked to this edge has already been explored")
        v = self._edge_to_vertex[e]
        v2 = self._edge_to_vertex[e2]
        geo = self._edge_to_geodesic[e]
        oc = self._explored[v]
        edges_to_flip = oc._correspondance[geo]
        tri = oc._triangulation.copy()
        
        for elt in edges_to_flip:
            tri.flip(elt)
        new_c = IsoDelaunayCell(tri)
        self._boundary.remove(v2)
        self._explored[v2] = new_c
        
        fp = self._graph.face_permutation(copy = False)
        #ep = self._graph.edge_permutation(copy = False)
        vp = self._graph.vertex_permutation(copy = False)
        new_ver = self._graph.num_vertices()
        new_edges = new_c.edges()
        nb_e = len(new_edges)
        for i in range(nb_e):
            if new_edges[i].unoriented() == geo.unoriented():
                index = i
        
        next_edge = vp[e2]
        pre_edge = e2
        self._edge_to_geodesic[e2] = geo
        new_boundary_face = True
        j = (index + 1) % nb_e
        cell_to_index = {}
        augment = False
        while next_edge != e2 or j != index: # we rotate around the vertex v2
            v_test = self._edge_to_vertex[2*(next_edge//2)+(1- next_edge%2)]
            c_test = self._explored[v_test]
            ng = self._edge_to_geodesic[2*(next_edge//2)+(1- next_edge%2)]
            cell_to_index[v_test] = True
            if new_c.adjacent(c_test) and ng.unoriented() == new_c._polygon.intersection(c_test._polygon).geodesic() and cell_to_index[v_test]:
                if new_edges[j].unoriented() == ng.unoriented():
                    self._edge_to_geodesic[next_edge] = ng
                    new_boundary_face = True
                    if augment:
                        new_ver += 1
                    pre_edge = next_edge
                    next_edge = vp[next_edge]
                    j = (j + 1) % nb_e
                    cell_to_index[v_test] = False
                else:
                    augment = True
                    if new_boundary_face:
                        self._edge_to_vertex.append(v2)
                        self._edge_to_vertex.append(new_ver)
                        self._edge_to_geodesic.append(new_edges[j])
                        self._edge_to_geodesic.append(None)
                        new_boundary_face = False
                        self._boundary.append(new_ver)
                        pre_edge = self._graph.add_edge(next_edge)
                        self._graph._check()
                    else:
                        self._edge_to_vertex.append(new_ver)
                        self._edge_to_vertex.append(v2)
                        self._edge_to_geodesic.append(None)
                        self._edge_to_geodesic.append(new_edges[j])
                        self._graph._realloc(self._graph._n+2)
                        self._graph.split_face(pre_edge,ner)
                        nel = fp[2*(nel//2)+(1- nel%2)]
                        pre_edge = vp[pre_edge]
                    j = (j + 1) % nb_e
            else:
                augment = True
                if new_boundary_face:
                    new_boundary_face = False
                    self._boundary.append(new_ver)
                    self._graph.move_dart(next_edge)
                    self._edge_to_vertex[next_edge] = new_ver
                    ner = next_edge # right next edge
                    nel = next_edge # left next edge
                    next_edge = vp[pre_edge]
                else:
                    self._graph.move_dart(next_edge, ner)
                    self._edge_to_vertex[next_edge] = new_ver
                    ner = next_edge
                    next_edge = vp[pre_edge]
        
    def is_explorable(self,e):
        r"""
        Test whether the edge ``e`` is explorable.

        EXAMPLES :

            sage: from flatsurf import translation_surfaces
            sage: S = translation_surfaces.arnoux_yoccoz(3)
            sage: t = FlatTriangulation(S.triangulate())
            sage: R = IsoDelaunayCell(t)
            sage: T = IsoDelaunayTessellation(R)
            sage: T.is_explorable(0)
            True
            sage: T.is_explorable(1)
            False
        """
        if self._edge_to_geodesic[e] == None:
            return False
        if e%2 == 0:
            e2 = e+1
        else : e2 = e-1
        if self._edge_to_geodesic[e2] != None:
            return False
        return True

    def explorable_edges(self):
        r"""
        Return a tuple containing the list of the label of explorable edges and the list of the corresponding geodesics.
        """
        res = []
        res2 = []
        for e in range(len(self._graph.vertex_permutation(copy=False))):
            if self.is_explorable(e):
                res.append(e)
                res2.append(self._edge_to_geodesic[e])
        return res,res2
            
        
    def plot(self, *args, **kwds):
        res = Graphics()
        for elt in self._explored.values():
            res = res + elt.plot(*args, **kwds)
        return res




