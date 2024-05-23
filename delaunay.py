"""

Iso-Delaunay

"""


from flatsurf import *
H = HyperbolicPlane(AA)
from surface_dynamics import *
from flatsurf.geometry.euclidean import ccw
# our first set of import. We may want to change it later


class FlatTriangulation:
    def __init__(self,M):
        if not M.is_triangulated():
            raise ValueError("The surface is not triangulated")
        self.surface = M
        edges = M.edges()
        m = len(edges)
        fp = [-1]*m
        ep = [2*(i//2)+(1- i%2) for i in range(m)]
        self.vectors = [] # the list of all the vectors of the edges : it contains the information of the holonomy and of the triangle (from the initial triangulation of M) it belongs
        labelling = {}
        for (a,b) in edges:
            if (a,b) not in labelling:
                (c,d) = M.opposite_edge(a,b)
                labelling[a,b]=len(labelling)
                self.vectors.append(M.tangent_vector(a,M.polygon(a).vertex(b),M.polygon(a).edge(b)))
                labelling[c,d]=len(labelling)
                self.vectors.append(M.tangent_vector(c,M.polygon(c).vertex(d),M.polygon(c).edge(d)))
        for lab in M.labels():
            for i in range(3):
                fp[labelling[lab, i]]=labelling[lab, (i+1)%3]

        # NOTE: on recent versions of surface_dynamics, ep is forced to be
        # 0 <-> 1, 2 <-> 3, etc and is not anymore an argument of the constructor
        try:
            self.triangulation = FatGraph(ep=ep, fp=fp) # the graph given by the triangulation
        except TypeError:
            self.triangulation = FatGraph(fp=fp) # the graph given by the triangulation

    def _check(self):
        self.triangulation._check()

        fp = self.triangulation.face_permutation(copy=False)
        n = len(self.vectors)
        for e1 in range(0, n, 2*n):
            v1 = self.vectors[e1].vector()
            if v1 != -self.vectors[e1 + 1].vector():
                raise ValueError('edge ({}, {}) have vectors with non-opposite holonomies'.format(e1, e1+1))

            e2 = fp[e1]
            v2 = self.vectors[e2].vector()
            e3 = fp[e2]
            v3 = self.vectors[e3].vector()
            if not (v1 + v2 + v3).is_zero() or ccw(v1, v2) <= 0 or ccw(v2, v3) <= 0 or ccw(v3, v1) <= 0:
                raise ValueError('({}, {}, {}) is not a valid triangle'.format(e1, e2, e3))

    def __repr__(self):
        return 'Flat Triangulation of a {}'.format(self.surface)

    def copy(self):
        nt = FlatTriangulation(self.surface)
        nt.vectors = self.vectors[:]
        nt.triangulation = self.triangulation.copy()
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
        if self.surface is not t.surface:
            raise ValueError('can not compare FlatTriangulation on different underlying surfaces')

        # find the first vector of self in other
        n = len(self.vectors)
        u0 = self.vectors[0]
        i = 0
        while i < n:
            v = other.vectors[i]
            if u0 == v:
                break
            i += 1
        if i == n:
            return (False, None) if certificate else False

        # does simultaneous DFS
        sfp = self.triangulation.face_permutation(copy=False)
        ofp = other.triangulation.face_permutation(copy=False)
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
                    if self.vectors[i] != other.vectors[ii]:
                        return (False, None) if certificate else False
                    relabelling[i] = ii
                    relabelling[i ^ 1] = ii ^ 1
                    todo.append(i ^ 1)
                elif relabelling[i] != ii:
                    return (False, None) if certificate else False
                i = sfp[i]
                ii = ofp[ii]

        return (True, relabelling) if certificate else True

    def is_flipable(self,e1):
        if e1%2 == 0:
            e2 = e1+1
        else : e2 = e1-1
        fp = self.triangulation.face_permutation(copy = False)
        f1 = fp[e1]
        f2 = fp[f1]
        g1 = fp[e2]
        g2 = fp[g1]
        return (ccw(self.vectors[f2].vector(),self.vectors[g1].vector()) > 0) and (ccw(self.vectors[g2].vector(),self.vectors[f1].vector()) > 0)

    def flip(self, e1):
        M = self.surface
        if not self.is_flipable(e1):
             raise ValueError("This edge is not flipable")
        # modification of the graph
        if e1%2 == 0:
            e2 = e1+1
        else : e2 = e1-1
        fp = self.triangulation.face_permutation(copy=True)
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
            self.triangulation = FatGraph(ep=ep, fp=fp)
        except TypeError:
            self.triangulation = FatGraph(fp=fp)

        # modification of the vectors
        hol2 = self.vectors[g1].vector()+self.vectors[f2].vector()
        hol1 = -hol2
        
        lab_pol = self.vectors[g2].polygon_label()
        lab_ver = self.vectors[g2].vertex()
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


        lab_pol = self.vectors[f2].polygon_label()
        lab_ver = self.vectors[f2].vertex()
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

        self.vectors[e1] = M.tangent_vector(pol1, M.polygon(pol1).vertex(ver1), hol1)
        self.vectors[e2] = M.tangent_vector(pol2, M.polygon(pol2).vertex(ver2), hol2)

    
    def geodesic(self, e1): #return the hyperbolic geometric associate to the edge with label e1 in the corresponding triangulation
        if not self.is_flipable(e1):
            return None
        if e1%2 == 0:
            e2 = e1+1
        else : e2 = e1-1
        fp = self.triangulation.face_permutation(copy = False)
        f1 = fp[e1]
        f2 = fp[f1]
        g1 = fp[e2]
        g2 = fp[g1]
        x1 = t.vectors[g1].vector()[0]
        y1 = t.vectors[g1].vector()[1]
        x2 = t.vectors[e1].vector()[0]
        y2 = t.vectors[e1].vector()[1]
        V = -t.vectors[f2].vector()
        x3 = V[0]
        y3 = V[1]
        a = x1*y2*y3*(y3-y2) + x2*y1*y3*(y1-y3) + x3*y1*y2*(y2-y1)
        b = x1*y1*(x2*y3-x3*y2) + x2*y2*(x3*y1-x1*y3) + x3*y3*(x1*y2-x2*y1)
        c = x1*x2*y3*(x1-x2) + x2*x3*y1*(x2-x3) + x1*x3*y2*(x3-x1)
        return H.geodesic(a,2*b,c=c,model = "half_plane")
    
    def plot(self, *args, **kwds):
        res = Graphics()
        for v in self.vectors[::2]:
            traj = v.straight_line_trajectory()
            traj.flow(100)
            res = res + traj.plot(*args, **kwds)
        return res





class IsoDelaunayCell:
    def __init__(self,t):
        self.triangulation = t.copy()
        Passage1 = True
        n = len(t.vectors)
        fp = t.triangulation.face_permutation(copy = False)
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
        self.polygon = Sol
        e = [s.geodesic() for s in Sol.edges()]
        cor_final = {}
        for g in e:
            cor_final[g] = cor[g]
        self.correspondance = cor_final

    def plot(self, *args, **kwds):
        return self.polygon.plot(*args, **kwds)





        