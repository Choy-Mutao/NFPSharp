namespace NFPSharp
{
    public struct TouchType
    {
        public int type { get; set; }
        public int A { get; set; }
        public int B { get; set; }
        public TouchType(int t, int a, int b)
        {
            type = t;
            A = a;
            B = b;
        }
    }

    public class NFPPoint
    {
        public double x, y;
        public bool marked;

        public NFPPoint() : this(0, 0) { }

        public NFPPoint(double _x, double _y)
        {
            x = _x;
            y = _y;
        }

        public static bool operator ==(NFPPoint p1, NFPPoint p2) => p1.x == p2.x && p1.y == p2.y;

        public static bool operator !=(NFPPoint p1, NFPPoint p2) => !(p1 == p2);
    }

    public class NFPVector
    {
        public double x, y;
        public NFPPoint start;
        public NFPPoint end;

        public NFPVector() : this(0, 0) { }
        public NFPVector(double _x, double _y, NFPPoint s, NFPPoint e)
        {
            x = _x; y = _y;
            start = s; end = e;
        }

        public NFPVector(double _x, double _y)
        {
            x = _x;
            y = _y;
            start = new NFPPoint(0, 0);
            end = new NFPPoint(x, y);
        }
        public NFPVector(NFPPoint s, NFPPoint e)
        {
            x = e.x - s.x;
            y = e.y - s.y;
            start = s; end = e;
        }

        public void SetLength(double D)
        {
            double unity = (x * x + y * y);
            if (unity == 0) throw new DivideByZeroException("This vector is zero");
            double d = Math.Sqrt(D * D / unity);
            x *= d;
            y *= d;
        }
    }

    public class NFPPolygon : List<NFPPoint>
    {
        private NFPPolygon() : base() { }

        public NFPPolygon(IEnumerable<double[]> double_array)
        {
            if (double_array.Any(item => item.Length != 2)) throw new ArgumentException("NFPPolygon constructor by double array is not valid");
            AddRange(double_array.Select(item => new NFPPoint(item[0], item[1])));
        }
        public NFPPolygon(IEnumerable<NFPPoint> points) : base(points) { }

        public double offsetx = 0;
        public double offsety = 0;
    }

    public static class MathUtil
    {
        public static bool AlmostEqual(double a, double b, double tolerance = 1e-9) => Math.Abs(a - b) < tolerance;

        /// <summary>
        /// </summary>
        /// <param name="p1"></param>
        /// <param name="p2"></param>
        /// <param name="distance"></param>
        /// <returns> true if points are within the given distances </returns>
        static bool WithinDistance(NFPPoint p1, NFPPoint p2, double distance)
        {
            var dx = p1.x - p2.x;
            var dy = p1.y - p2.y;
            return ((dx * dx + dy * dy) < distance * distance);
        }

        static double DegreesToRadians(double angle) => angle * (Math.PI / 180);

        static double RadiansToDegrees(double angle) => angle * (180 / Math.PI);

        /// <summary>
        /// normalize vector into a unit vector
        /// </summary>
        /// <param name="v"></param>
        /// <returns></returns>
        public static NFPVector NormalizeVector(NFPVector v)
        {
            if (AlmostEqual(v.x * v.x + v.y * v.y, 1))
                return v; // given vector was already a unit vector
            var len = Math.Sqrt(v.x * v.x + v.y * v.y);
            var inverse = 1 / len;

            double x = v.x * inverse, y = v.y * inverse;
            return new NFPVector(x, y);
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="A"></param>
        /// <param name="B"></param>
        /// <param name="p"></param>
        /// <returns>
        ///  true if p lies on the line segment defined by AB, but <b>not at any endpoints</b>
        /// </returns>
        public static bool OnSegment(NFPPoint A, NFPPoint B, NFPPoint p, double tolerance)
        {
            double lengthSqr1 = (p.x - A.x) * (p.x - A.x) + (p.y - A.y) * (p.y - A.y);
            double lengthSqr2 = (p.x - B.x) * (p.x - B.x) + (p.y - B.y) * (p.y - B.y);
            if (AlmostEqual(lengthSqr1, 0, tolerance) || AlmostEqual(lengthSqr2, 0, tolerance)) return false;

            double lengthSqr_ref = (A.x - B.x) * (A.x - B.x) + (A.y - B.y) * (A.y - B.y);
            double l1 = Math.Sqrt(lengthSqr1);
            double l2 = Math.Sqrt(lengthSqr2);
            double l_ref = Math.Sqrt(lengthSqr_ref);

            if (!AlmostEqual(l1 + l2, l_ref, tolerance)) return false;
            return true;
        }
    }

    public class NFPUtils
    {
        /// <summary>
        /// given a static polygon A and a movable polygon B, compute a no fit polygon by orbiting B about A
        /// if the inside flag is set, B is orbited inside of A rather than outside
        /// if the searchEdges flag is set, all edges of A are explored for NFPs - multiple 
        /// </summary>
        /// <param name="A"></param>
        /// <param name="B"></param>
        /// <param name="inside"></param>
        /// <param name="searchEdges"></param>
        /// <returns></returns>
        public static List<NFPPolygon>? NoFitPolygon(NFPPolygon A, NFPPolygon B, bool inside, bool searchEdges, double tolerance = 1e-9)
        {
            if (A.Count < 3 || B.Count < 3) return null;

            A.offsetx = 0;
            A.offsety = 0;

            int minAindex = 0; double minA = A[minAindex].y;
            int maxBindex = 0; double maxB = B[maxBindex].y;

            for (int i = 1; i < A.Count; i++)
            {
                A[i].marked = false;
                if (A[i].y < minA)
                {
                    minA = A[i].y;
                    minAindex = i;
                }
            }

            for (int i = 1; i < B.Count; i++)
            {
                B[i].marked = false;
                if (B[i].y > maxB)
                {
                    maxB = B[i].y;
                    maxBindex = i;
                }
            }

            NFPPoint? startpoint;
            if (!inside)
            {
                // Shift B such that the bottom-most point of B is at the top-most point of A.
                // This guarantees an initial placement with no intersections
                startpoint = new NFPPoint(A[minAindex].x - B[maxBindex].x,
                                            A[minAindex].y - B[maxBindex].y);
            }
            else
            {
                // no reliable heuristic for inside
                startpoint = SearchStartPoint(A, B, true);
            }

            List<NFPPolygon> NFPlist = [];
            while (startpoint is not null)
            {
                B.offsetx = startpoint.x;
                B.offsety = startpoint.y;

                NFPVector? prevvector = null; // keep track of previous vector
                NFPPolygon NFP = new([[B[0].x + B.offsetx, B[0].y + B.offsety]]);

                double referencex = B[0].x + B.offsetx;
                double referencey = B[0].y + B.offsety;
                double startx = referencex;
                double starty = referencey;

                double counter = 0;
                while (counter < 10 * (A.Count + B.Count)) // sanity check, prevent inifinite loop
                {
                    var touching = FindTouching(A, B, tolerance);
                    var vectors = GenerateTranslationVectors(touching, ref A, ref B).ToArray();

                    // todo: there should be a faster way to reject vectors that will cause immediate intersection. For now just check them all

                    NFPVector? translate = null;
                    double maxd = 0;
                    for (int i = 0; i < vectors.Length; i++)
                    {
                        if (vectors[i].x == 0 && vectors[i].y == 0)
                            continue;

                        // if this vector points us back to where we came from, ignore it.
                        // ie cross product = 0, dot product < 0
                        if (prevvector != null && vectors[i].y * prevvector.y + vectors[i].x * prevvector.x < 0)
                        {
                            // compare magnitude with unity vectors
                            double vectorlength = Math.Sqrt(vectors[i].x * vectors[i].x + vectors[i].y * vectors[i].y);
                            NFPVector unitv = new NFPVector(vectors[i].x / vectorlength, vectors[i].y / vectorlength);

                            double prevlength = Math.Sqrt(prevvector.x * prevvector.x + prevvector.y * prevvector.y);
                            NFPVector prevunit = new NFPVector(prevvector.x / prevlength, prevvector.y / prevlength);

                            // we need to scale down to unit vectors to normalize vector length. Could also just do a tan here
                            if (Math.Abs(unitv.y * prevunit.x - unitv.x * prevunit.y) < 0.0001)
                            {
                                continue;
                            }
                        }

                        var d = PolygonSlideDistance(A, B, vectors[i], true);
                        var vecd2 = vectors[i].x * vectors[i].x + vectors[i].y * vectors[i].y;

                        if (d == null || d.Value * d.Value > vecd2)
                        {
                            var vecd = Math.Sqrt(vectors[i].x * vectors[i].x + vectors[i].y * vectors[i].y);
                            d = vecd;
                        }

                        if (d != null && d > maxd)
                        {
                            maxd = d.Value;
                            translate = vectors[i];
                        }
                    }

                    if (translate == null || MathUtil.AlmostEqual(maxd, 0))
                    {
                        // didn't close the loop, something went wrong here
                        NFP = null;
                        break;
                    }

                    translate.start.marked = true;
                    translate.end.marked = true;

                    prevvector = translate;

                    // trim
                    double vlength2 = translate.x * translate.x + translate.y * translate.y;
                    if (maxd * maxd < vlength2 && MathUtil.AlmostEqual(maxd * maxd, vlength2))
                    {
                        var scale = Math.Sqrt((maxd * maxd) / vlength2);
                        translate.x *= scale;
                        translate.y *= scale;
                    }

                    referencex += translate.x;
                    referencey += translate.y;

                    if (MathUtil.AlmostEqual(referencex, startx) && MathUtil.AlmostEqual(referencey, starty))
                    {
                        // we've made a full loop
                        break;
                    }

                    // if A and B start on a touching horizontal line, the end point may not be the start point
                    bool looped = false;
                    if (NFP.Count > 0)
                    {
                        for (int i = 0; i < NFP.Count - 1; i++)
                        {
                            if (MathUtil.AlmostEqual(referencex, NFP[i].x) && MathUtil.AlmostEqual(referencey, NFP[i].y))
                            {
                                looped = true;
                            }
                        }
                    }

                    if (looped)
                    {
                        // we've made a full loop
                        break;
                    }

                    NFP.Add(new NFPPoint(referencex, referencey));

                    B.offsetx += translate.x;
                    B.offsety += translate.y;

                    counter++;
                }

                if (NFP != null && NFP.Count > 0) NFPlist.Add(NFP);

                if (!searchEdges) break; // only get outer NFP or first inner NFP;

                startpoint = SearchStartPoint(A, B, inside);
            }

            return NFPlist;
        }

        /// <summary>
        /// Find touching vertices/edges
        /// </summary>
        /// <param name="A"></param>
        /// <param name="B"></param>
        /// <param name="tolerance"></param>
        /// <returns></returns>
        private static IEnumerable<TouchType> FindTouching(NFPPolygon A, NFPPolygon B, double tolerance)
        {
            List<TouchType> touching = new List<TouchType>();
            for (int i = 0; i < A.Count; i++)
            {
                int nexti = (i == A.Count - 1) ? 0 : i + 1;
                for (int j = 0; j < B.Count; j++)
                {
                    int nextj = (j == B.Count - 1) ? 0 : j + 1;
                    if (MathUtil.AlmostEqual(A[i].x, B[j].x + B.offsetx) && MathUtil.AlmostEqual(A[i].y, B[j].y + B.offsety))
                    {
                        touching.Add(new TouchType(0, i, j));
                    }
                    else if (MathUtil.OnSegment(A[i], A[nexti], new NFPPoint(B[j].x + B.offsetx, B[j].y + B.offsety), tolerance))
                    {
                        touching.Add(new TouchType(1, nexti, j));
                    }
                    else if (true)
                    {
                        touching.Add(new TouchType(2, i, nextj));
                    }
                }
            }
            return touching;
        }

        private static IEnumerable<NFPVector> GenerateTranslationVectors(IEnumerable<TouchType> touching, ref NFPPolygon A, ref NFPPolygon B)
        {
            List<NFPVector> vectors = new List<NFPVector>();
            for (int i = 0; i < touching.Count(); i++)
            {
                var vertexA = A[touching.ElementAt(i).A];
                vertexA.marked = true;

                // adjacent A vertices
                var prevAindex = touching.ElementAt(i).A - 1;
                var nextAindex = touching.ElementAt(i).A + 1;

                prevAindex = (prevAindex < 0) ? A.Count - 1 : prevAindex; // loop
                nextAindex = (nextAindex >= A.Count) ? 0 : nextAindex; // loop

                var prevA = A[prevAindex];
                var nextA = A[nextAindex];

                // adjacent B vertices
                var vertexB = B[touching.ElementAt(i).B];

                var prevBindex = touching.ElementAt(i).B - 1;
                var nextBindex = touching.ElementAt(i).B + 1;

                prevBindex = (prevBindex < 0) ? B.Count - 1 : prevBindex; // loop
                nextBindex = (nextBindex >= B.Count) ? 0 : nextBindex; // loop

                var prevB = B[prevBindex];
                var nextB = B[nextBindex];

                if (touching.ElementAt(i).type == 0)
                {
                    NFPVector vA1 = new NFPVector(vertexA, prevA), vA2 = new NFPVector(vertexA, nextA);
                    // B vectors need to be inverted
                    NFPVector vB1 = new NFPVector(prevB, vertexB), vB2 = new NFPVector(nextB, vertexB);


                    vectors.Add(vA1);
                    vectors.Add(vA2);
                    vectors.Add(vB1);
                    vectors.Add(vB2);
                }
                else if (touching.ElementAt(i).type == 1)
                {
                    vectors.Add(new NFPVector(
                        vertexA.x - (vertexB.x + B.offsetx),
                        vertexA.y - (vertexB.y + B.offsety),
                        prevA,
                        vertexA
                        ));

                    vectors.Add(new NFPVector(
                        prevA.x - (vertexB.x + B.offsetx),
                        prevA.y - (vertexB.y + B.offsety),
                        vertexA,
                        prevA
                        ));
                }
                else if (touching.ElementAt(i).type == 2)
                {
                    vectors.Add(new NFPVector(
                        vertexA.x - (vertexB.x + B.offsetx),
                        vertexA.y - (vertexB.y + B.offsety),
                        prevB,
                        vertexB
                        ));

                    vectors.Add(new NFPVector(
                        vertexA.x - (prevB.x + B.offsetx),
                        vertexA.y - (prevB.y + B.offsety),
                        vertexB,
                        prevB
                        ));
                }
            }
            return vectors;
        }

        /// <summary>
        /// Searches for an arrangement of A and B such that they do not overlap
        /// if an NFP is given, only search for start points that have not already been traversed in the given NFP
        /// </summary>
        /// <param name="A"></param>
        /// <param name="B"></param>
        /// <param name="inside"></param>
        /// <returns></returns>
        /// <exception cref="NotImplementedException"></exception>
        public static NFPPoint? SearchStartPoint(in NFPPolygon _A, in NFPPolygon _B, bool inside, List<NFPPolygon>? NFP = null, double tolerance = 1e-9)
        {
            // Clone Arrays
            NFPPolygon A = new(_A);
            NFPPolygon B = new(_B);

            // close the loop for polygons
            if (A.First() != A.Last()) A.Add(A.First());
            if (B.First() != B.Last()) B.Add(B.First());

            for (int i = 0; i < A.Count - 1; i++)
            {
                if (!A[i].marked)
                {
                    A[i].marked = true;
                    for (int j = 0; j < B.Count; j++)
                    {
                        B.offsetx = A[i].x - B[j].x;
                        B.offsety = A[i].y - B[j].y;

                        bool? Binside = null;
                        for (int k = 0; k < B.Count; k++)
                        {
                            var inpoly = PointInPolygon(new NFPPoint(B[k].x + B.offsetx, B[k].y + B.offsety), A, tolerance);
                            if (inpoly != null)
                            {
                                Binside = inpoly;
                                break;
                            }
                        }

                        if (Binside == null) return null; // A and B are the same
                        var startPoint = new NFPPoint(B.offsetx, B.offsety);
                        if (((Binside.Value == true && inside) || (Binside.Value == false && !inside)) && !Intersect(A, B, tolerance) && !InNfp(startPoint, NFP, tolerance)) return startPoint;

                        // slide B along vector
                        var vx = A[i + 1].x - A[i].x;
                        var vy = A[i + 1].y - A[i].y;

                        double? d1 = polygonProjectionDistance(A, B, new NFPPoint(vx, vy));
                        double? d2 = polygonProjectionDistance(A, B, new NFPPoint(-vx, -vy));

                        double? d = null;

                        //todo: clean this up
                        if (d1 is null && d2 is null) { /* nothing */ }
                        else if (d1 is null) { d = d2; }
                        else if (d2 is null) { d = d1; }
                        else { d = Math.Min(d1.Value, d2.Value); }

                        // only slide until no longer negative
                        // todo: clean this up
                        if (d is not null && MathUtil.AlmostEqual(d.Value, 0, tolerance) && d.Value > 0)
                        {
                            // nothing
                        }
                        else continue;

                        var vd2 = vx * vx + vy * vy;

                        if((d.Value * d.Value) < vd2 && !MathUtil.AlmostEqual(d.Value * d.Value, vd2))
                        {
                            var vd = Math.Sqrt(vx * vx + vy * vy);
                            vx *= d.Value / vd;
                            vy *= d.Value / vd;
                        }

                        B.offsetx += vx;
                        B.offsety += vy;

                        for(int k = 0; k < B.Count; k++)
                        {
                            var inpoly = PointInPolygon(new NFPPoint(B[k].x + B.offsetx, B[k].y + B.offsety), A, tolerance);
                            if(inpoly != null)
                            {
                                Binside = inpoly;
                                break;
                            }
                        }
                        startPoint = new NFPPoint(B.offsetx, B.offsety);
                        if (((Binside.Value == true && inside) || (Binside.Value == false && !inside)) && !Intersect(A, B, tolerance) && !InNfp(startPoint, NFP, tolerance))
                        {
                            return startPoint;
                        }
                    }
                }
            }

            return null;
        }

        /// <summary>
        /// return true if point already exists in the given nfps
        /// </summary>
        /// <param name="p"></param>
        /// <param name="nfp"></param>
        /// <returns></returns>
        private static bool InNfp(NFPPoint p, List<NFPPolygon> nfp, double tolerance)
        {
            if (nfp == null || nfp.Count == 0)
                return false;

            for (var i = 0; i < nfp.Count; i++)
                for (var j = 0; j < nfp[i].Count; j++)
                    if (MathUtil.AlmostEqual(p.x, nfp[i][j].x, tolerance)
                        && MathUtil.AlmostEqual(p.y, nfp[i][j].y, tolerance))
                        return true;

            return false;
        }

        /// <summary>
        /// TODO: swap this for a more efficient sweep-line implementation
        /// return Edges: if SetByRowElements, return all edges on A that have intersections
        /// </summary>
        /// <param name="A"></param>
        /// <param name="B"></param>
        /// <returns></returns>
        private static bool Intersect(in NFPPolygon A, in NFPPolygon B, double tolerance)
        {
            var Aoffsetx = A.offsetx;
            var Aoffsety = A.offsety;

            var Boffsetx = B.offsetx;
            var Boffsety = B.offsety;

            //A = A.slice(0);
            //B = B.slice(0);

            for (int i = 0; i < A.Count - 1; i++)
            {
                for (int j = 0; j < B.Count - 1; j++)
                {
                    NFPPoint a1 = new NFPPoint(A[i].x + Aoffsetx, A[i].y + Aoffsety),
                        a2 = new NFPPoint(A[i + 1].x + Aoffsetx, A[i + 1].y + Aoffsety),
                        b1 = new NFPPoint(B[j].x + Boffsetx, B[j].y + Boffsety),
                        b2 = new NFPPoint(B[j + 1].x + Boffsetx, B[j + 1].y + Boffsety);

                    var prevbindex = (j == 0) ? B.Count - 1 : j - 1;
                    var prevaindex = (i == 0) ? A.Count - 1 : i - 1;
                    var nextbindex = (j + 1 == B.Count - 1) ? 0 : j + 2;
                    var nextaindex = (i + 1 == A.Count - 1) ? 0 : i + 2;

                    // go even further back if we happen to hit on a loop end point
                    if (B[prevbindex] == B[j] || (MathUtil.AlmostEqual(B[prevbindex].x, B[j].x, tolerance) && MathUtil.AlmostEqual(B[prevbindex].y, B[j].y, tolerance)))
                    {
                        prevbindex = (prevbindex == 0) ? B.Count - 1 : prevbindex - 1;
                    }

                    if (A[prevaindex] == A[i] || (MathUtil.AlmostEqual(A[prevaindex].x, A[i].x, tolerance) && MathUtil.AlmostEqual(A[prevaindex].y, A[i].y, tolerance)))
                    {
                        prevaindex = (prevaindex == 0) ? A.Count - 1 : prevaindex - 1;
                    }

                    // go even further forward if we happen to hit on a loop end point
                    if (B[nextbindex] == B[j + 1] || (MathUtil.AlmostEqual(B[nextbindex].x, B[j + 1].x, tolerance) && MathUtil.AlmostEqual(B[nextbindex].y, B[j + 1].y, tolerance)))
                    {
                        nextbindex = (nextbindex == B.Count - 1) ? 0 : nextbindex + 1;
                    }

                    if (A[nextaindex] == A[i + 1] || (MathUtil.AlmostEqual(A[nextaindex].x, A[i + 1].x, tolerance) && MathUtil.AlmostEqual(A[nextaindex].y, A[i + 1].y, tolerance)))
                    {
                        nextaindex = (nextaindex == A.Count - 1) ? 0 : nextaindex + 1;
                    }

                    NFPPoint a0 = new NFPPoint(A[prevaindex].x + Aoffsetx, A[prevaindex].y + Aoffsety),
                        b0 = new NFPPoint(B[prevbindex].x + Boffsetx, B[prevbindex].y + Boffsety),
                        a3 = new NFPPoint(A[nextaindex].x + Aoffsetx, A[nextaindex].y + Aoffsety),
                        b3 = new NFPPoint(B[nextbindex].x + Boffsetx, B[nextbindex].y + Boffsety);


                    if (MathUtil.OnSegment(a1, a2, b1, tolerance) || (MathUtil.AlmostEqual(a1.x, b1.x, tolerance) && MathUtil.AlmostEqual(a1.y, b1.y, tolerance)))
                    {
                        // if a point is on a segment, it could intersect or it could not. Check via the neighboring points
                        var b0in = PointInPolygon(b0, A, tolerance);
                        var b2in = PointInPolygon(b2, A, tolerance);

                        // todo: b0in, b2in maybe null
                        if ((b0in == true && b2in == false) || (b0in == false && b2in == true))
                        {
                            return true;
                        }
                        else
                        {
                            continue;
                        }
                    }

                    if (MathUtil.OnSegment(a1, a2, b2, tolerance) || (MathUtil.AlmostEqual(a2.x, b2.x, tolerance) && MathUtil.AlmostEqual(a2.y, b2.y, tolerance)))
                    {
                        // if a point is on a segment, it could intersect or it could not. Check via the neighboring points
                        var b1in = PointInPolygon(b1, A, tolerance);
                        var b3in = PointInPolygon(b3, A, tolerance);

                        // todo: b0in, b2in maybe null
                        if ((b1in == true && b3in == false) || (b1in == false && b3in == true))
                        {
                            return true;
                        }
                        else
                        {
                            continue;
                        }
                    }

                    if (MathUtil.OnSegment(b1, b2, a1, tolerance) || (MathUtil.AlmostEqual(a1.x, b2.x, tolerance) && MathUtil.AlmostEqual(a1.y, b2.y, tolerance)))
                    {
                        // if a point is on a segment, it could intersect or it could not. Check via the neighboring points
                        var a0in = PointInPolygon(a0, B, tolerance);
                        var a2in = PointInPolygon(a2, B, tolerance);

                        // todo: b0in, b2in maybe null
                        if ((a0in == true && a2in == false) || (a0in == false && a2in == true))
                        {
                            return true;
                        }
                        else
                        {
                            continue;
                        }
                    }

                    if (MathUtil.OnSegment(b1, b2, a2, tolerance) || (MathUtil.AlmostEqual(a2.x, b1.x, tolerance) && MathUtil.AlmostEqual(a2.y, b1.y, tolerance)))
                    {
                        // if a point is on a segment, it could intersect or it could not. Check via the neighboring points
                        var a1in = PointInPolygon(a1, B, tolerance);
                        var a3in = PointInPolygon(a3, B, tolerance);

                        if ((a1in == true && a3in == false) || (a1in == false && a3in == true))
                        {
                            return true;
                        }
                        else
                        {
                            continue;
                        }
                    }
                    var p = LineIntersect(b1, b2, a1, a2);

                    if (p is not null)
                        return true;
                }
            }
            return false;
        }

        private static double? polygonProjectionDistance(NFPPolygon A, NFPPolygon B, NFPPoint p)
        {
            throw new NotImplementedException();
        }

        /// <summary>
        /// returns the intersection of AB and EF
        /// or null if there are no intersections or other numerical error
        /// if the infinite flag is set, AE and EF describe infinite lines without 
        /// endpoints, they are finite line segments otherwise
        /// </summary>
        /// <param name="A"></param>
        /// <param name="B"></param>
        /// <param name="E"></param>
        /// <param name="F"></param>
        /// <param name="infinite"></param>
        /// <param name="tolerance"></param>
        /// <returns></returns>
        public static NFPPoint? LineIntersect(NFPPoint A, NFPPoint B, NFPPoint E, NFPPoint F, bool infinite = false, double tolerance = 1e-9)
        {
            double a1, a2, b1, b2, c1, c2, x, y;

            a1 = B.y - A.y;
            b1 = A.x - B.x;
            c1 = B.x * A.y - A.x * B.y;
            a2 = F.y - E.y;
            b2 = E.x - F.x;
            c2 = F.x * E.y - E.x * F.y;

            var denom = a1 * b2 - a2 * b1;

            x = (b1 * c2 - b2 * c1) / denom;
            y = (a2 * c1 - a1 * c2) / denom;

            if (double.IsFinite(x) || double.IsFinite(y)) return null;

            if (!infinite)
            {
                // coincident points do not count as intersecting
                if (Math.Abs(A.x - B.x) > tolerance && ((A.x < B.x) ? x < A.x || x > B.x : x > A.x || x < B.x)) return null;
                if (Math.Abs(A.y - B.y) > tolerance && ((A.y < B.y) ? y < A.y || y > B.y : y > A.y || y < B.y)) return null;

                if (Math.Abs(E.x - F.x) > tolerance && ((E.x < F.x) ? x < E.x || x > F.x : x > E.x || x < F.x)) return null;
                if (Math.Abs(E.y - F.y) > tolerance && ((E.y < F.y) ? y < E.y || y > F.y : y > E.y || y < F.y)) return null;
            }

            return new NFPPoint(x, y);
        }

        private static bool? PointInPolygon(NFPPoint point, NFPPolygon polygon, double tolerance)
        {
            if (polygon.Count < 3)
            {
                return null;
            }

            var inside = false;
            var offsetx = polygon.offsetx;
            var offsety = polygon.offsety;

            for (int i = 0, j = polygon.Count - 1; i < polygon.Count; j = i++)
            {
                var xi = polygon[i].x + offsetx;
                var yi = polygon[i].y + offsety;
                var xj = polygon[j].x + offsetx;
                var yj = polygon[j].y + offsety;

                if (MathUtil.AlmostEqual(xi, point.x, tolerance) && MathUtil.AlmostEqual(yi, point.y, tolerance))
                {
                    return null; // no result
                }

                if (MathUtil.OnSegment(new NFPPoint(xi, yi), new NFPPoint(xj, yj), point, tolerance))
                {
                    return false; // exactly on the segment
                }

                if (MathUtil.AlmostEqual(xi, xj, tolerance) && MathUtil.AlmostEqual(yi, yj, tolerance))
                { // ignore very small lines
                    continue;
                }

                var intersect = ((yi > point.y) != (yj > point.y)) && (point.x < (xj - xi) * (point.y - yi) / (yj - yi) + xi);
                if (intersect) inside = !inside;
            }

            return inside;
        }

        public static double? PolygonSlideDistance(in NFPPolygon A, in NFPPolygon B, NFPVector direction, bool ignoreNegative = false, double tolerance = 1e-9)
        {
            double Aoffsetx, Aoffsety, Boffsetx, Boffsety;

            Aoffsetx = A.offsetx;
            Aoffsety = A.offsety;

            Boffsetx = B.offsetx;
            Boffsety = B.offsety;

            if (A.First() != A.Last()) A.Add(A.First());
            if (B.First() != B.Last()) B.Add(B.First());

            var edgeA = A;
            var edgeB = B;

            double? distance = null, d;

            var dir = MathUtil.NormalizeVector(direction);

            NFPPoint A1, A2, B1, B2;
            for (int i = 0; i < edgeB.Count - 1; i++)
            {
                for (int j = 0; j < edgeA.Count - 1; j++)
                {
                    A1 = new NFPPoint(edgeA[j].x + Aoffsetx, edgeA[j].y + Aoffsety);
                    A2 = new NFPPoint(edgeA[j + 1].x + A.offsetx, edgeA[j + 1].y + Aoffsety);

                    B1 = new NFPPoint(edgeB[i].x + Boffsetx, edgeB[i].y + Boffsety);
                    B2 = new NFPPoint(edgeB[i + 1].x + Boffsetx, edgeB[i + 1].y + Boffsety);

                    // ignore any one of lines if is too small
                    if ((MathUtil.AlmostEqual(A1.x, A2.x, tolerance) && MathUtil.AlmostEqual(A1.y, A2.y, tolerance)) || (MathUtil.AlmostEqual(B1.x, B2.x, tolerance) && MathUtil.AlmostEqual(B1.y, B2.y, tolerance))) continue;

                    d = SegmentDistance(A1, A2, B1, B2, dir, tolerance);

                    if ((d != null) && (distance == null || d.Value < distance))
                    {
                        if (!ignoreNegative || d.Value > 0 || MathUtil.AlmostEqual(d.Value, 0))
                            distance = d;
                    }
                }
            }

            return distance;
        }

        /// <summary>
        /// Calculate the distance between setment(A,B) and segment(E,F) in a certain direction
        /// </summary>
        /// <param name="A"></param>
        /// <param name="B"></param>
        /// <param name="E"></param>
        /// <param name="F"></param>
        /// <param name="direction"></param>
        /// <param name="tolerance"></param>
        /// <returns></returns>
        /// <exception cref="DivideByZeroException"></exception>
        public static double? SegmentDistance(NFPPoint A, NFPPoint B, NFPPoint E, NFPPoint F, NFPVector direction, double tolerance = 1e-9)
        {
            NFPVector normal = new NFPVector(direction.y, -direction.x);
            NFPVector reverse = new NFPVector(-direction.x, -direction.y);

            double dotA = A.x * normal.x + A.y * normal.y; // projection distance on noraml
            double dotB = B.x * normal.x + B.y * normal.y; // projection distance on noraml
            double dotE = E.x * normal.x + E.y * normal.y; // projection distance on noraml
            double dotF = F.x * normal.x + F.y * normal.y; // projection distance on noraml

            double crossA = A.x * direction.x + A.y * direction.y;
            double crossB = B.x * direction.x + B.y * direction.y;
            double crossE = E.x * direction.x + E.y * direction.y;
            double crossF = F.x * direction.x + F.y * direction.y;

            var ABmin = Math.Min(dotA, dotB);
            var ABmax = Math.Max(dotA, dotB);

            var EFmin = Math.Min(dotE, dotF);
            var EFmax = Math.Max(dotE, dotF);

            // segments that will merely touch at one point
            if (MathUtil.AlmostEqual(ABmax, EFmin, tolerance) || MathUtil.AlmostEqual(ABmin, EFmax, tolerance))
                return null;

            // segments miss eachother completely
            if (ABmax < EFmin || ABmin > EFmax) return null;

            int overlap;
            if ((ABmax > EFmax && ABmin < EFmin) || (EFmax > ABmax && EFmin < ABmin)) overlap = 1;
            else
            {
                var minMax = Math.Min(ABmax, EFmax);
                var maxMin = Math.Max(ABmin, EFmin);

                var maxMax = Math.Max(ABmax, EFmax);
                var minMin = Math.Min(ABmin, EFmin);

                if (maxMax - minMin == 0) throw new DivideByZeroException("SegmentDistance: maxMax - minMin == 0");
                overlap = (int)((minMax - maxMin) / (maxMax - minMin));
            }

            var crossABE = (E.y - A.y) * (B.x - A.x) - (E.x - A.x) * (B.y - A.y);
            var crossABF = (F.y - A.y) * (B.x - A.x) - (F.x - A.x) * (B.y - A.y);

            // lines are colinear
            if (MathUtil.AlmostEqual(crossABE, 0, tolerance) && MathUtil.AlmostEqual(crossABF, 0, tolerance))
            {

                NFPVector ABnorm = new NFPVector(B.y - A.y, A.x - B.x);
                NFPVector EFnorm = new NFPVector(F.y - E.y, E.x - F.x);

                var ABnormlength = Math.Sqrt(ABnorm.x * ABnorm.x + ABnorm.y * ABnorm.y);
                ABnorm.x /= ABnormlength;
                ABnorm.y /= ABnormlength;

                var EFnormlength = Math.Sqrt(EFnorm.x * EFnorm.x + EFnorm.y * EFnorm.y);
                EFnorm.x /= EFnormlength;
                EFnorm.y /= EFnormlength;

                // segment normals must point in opposite directions
                if (Math.Abs(ABnorm.y * EFnorm.x - ABnorm.x * EFnorm.y) < tolerance && ABnorm.y * EFnorm.y + ABnorm.x * EFnorm.x < 0)
                {
                    // normal of AB segment must point in same direction as given direction vector
                    var normdot = ABnorm.y * direction.y + ABnorm.x * direction.x;
                    // the segments merely slide along eachother
                    if (MathUtil.AlmostEqual(normdot, 0, tolerance))
                    {
                        return null;
                    }
                    if (normdot < 0)
                    {
                        return 0;
                    }
                }
                return null;
            }

            List<double> distances = [];

            // coincident points
            if (MathUtil.AlmostEqual(dotA, dotE)) distances.Add(crossA - crossE);
            else if (MathUtil.AlmostEqual(dotA, dotF)) distances.Add(crossA - crossF);
            else if (dotA > EFmin && dotA < EFmax)
            {
                var d = PointDistance(A, E, F, reverse);
                if (d != null && MathUtil.AlmostEqual(d.Value, 0))
                { //  A currently touches EF, but AB is moving away from EF
                    var dB = PointDistance(B, E, F, reverse, true);
                    if (dB < 0 || MathUtil.AlmostEqual(dB!.Value * overlap, 0))
                    {
                        d = null;
                    }
                }
                if (d != null)
                {
                    distances.Add(d.Value);
                }
            }

            if (MathUtil.AlmostEqual(dotB, dotE))
            {
                distances.Add(crossB - crossE);
            }
            else if (MathUtil.AlmostEqual(dotB, dotF))
            {
                distances.Add(crossB - crossF);
            }
            else if (dotB > EFmin && dotB < EFmax)
            {
                var d = PointDistance(B, E, F, reverse);

                if (d != null && MathUtil.AlmostEqual(d.Value, 0))
                { // crossA>crossB A currently touches EF, but AB is moving away from EF
                    var dA = PointDistance(A, E, F, reverse, true);
                    if (dA < 0 || MathUtil.AlmostEqual(dA!.Value * overlap, 0))
                    {
                        d = null;
                    }
                }
                if (d != null)
                {
                    distances.Add(d.Value);
                }
            }

            if (dotE > ABmin && dotE < ABmax)
            {
                var d = PointDistance(E, A, B, direction);
                if (d != null && MathUtil.AlmostEqual(d.Value, 0))
                { // crossF<crossE A currently touches EF, but AB is moving away from EF
                    var dF = PointDistance(F, A, B, direction, true);
                    if (dF < 0 || MathUtil.AlmostEqual(dF!.Value * overlap, 0))
                    {
                        d = null;
                    }
                }
                if (d != null)
                {
                    distances.Add(d.Value);
                }
            }

            if (dotF > ABmin && dotF < ABmax)
            {
                var d = PointDistance(F, A, B, direction);
                if (d != null && MathUtil.AlmostEqual(d.Value, 0))
                {
                    // && crossE<crossF A currently touches EF, but AB is moving away from EF
                    var dE = PointDistance(E, A, B, direction, true);
                    if (dE < 0 || MathUtil.AlmostEqual(dE!.Value * overlap, 0))
                    {
                        d = null;
                    }
                }
                if (d != null)
                {
                    distances.Add(d.Value);
                }
            }

            if (distances.Count == 0) return null;

            return distances.Min();
        }

        /// <summary>
        /// Find the distance from a point in the normal direction to a straight line
        /// If return is negative, it means in the opposite direction
        /// </summary>
        /// <param name="p"></param>
        /// <param name="s1"></param>
        /// <param name="s2"></param>
        /// <param name="normal"></param>
        /// <param name="infinite"></param>
        /// <returns></returns>
        public static double? PointDistance(NFPPoint p, NFPPoint s1, NFPPoint s2, NFPVector normal, bool infinite = false)
        {
            normal = MathUtil.NormalizeVector(normal);
            NFPVector dir = new NFPVector(normal.y, -normal.x);

            double pdot = p.x * dir.x + p.y * dir.y;
            double s1dot = s1.x * dir.x + s1.y * dir.y;
            double s2dot = s2.x * dir.x + s2.y * dir.y;

            double pdotnorm = p.x * normal.x + p.y * normal.y;
            double s1dotnorm = s1.x * normal.x + s1.y * normal.y;
            double s2dotnorm = s2.x * normal.x + s2.y * normal.y;

            if (!infinite)
            {
                if (((pdot < s1dot || MathUtil.AlmostEqual(pdot, s1dot)) && (pdot < s2dot || MathUtil.AlmostEqual(pdot, s2dot))) || ((pdot > s1dot || MathUtil.AlmostEqual(pdot, s1dot)) && (pdot > s2dot || MathUtil.AlmostEqual(pdot, s2dot))))
                {
                    return null; // dot doesn't collide with segment, or lies directly on the vertex
                }
                if ((MathUtil.AlmostEqual(pdot, s1dot) && MathUtil.AlmostEqual(pdot, s2dot)) && (pdotnorm > s1dotnorm && pdotnorm > s2dotnorm))
                {
                    return Math.Min(pdotnorm - s1dotnorm, pdotnorm - s2dotnorm);
                }
                if ((MathUtil.AlmostEqual(pdot, s1dot) && MathUtil.AlmostEqual(pdot, s2dot)) && (pdotnorm < s1dotnorm && pdotnorm < s2dotnorm))
                {
                    return -Math.Min(s1dotnorm - pdotnorm, s2dotnorm - pdotnorm);
                }
            }

            return -(pdotnorm - s1dotnorm + (s1dotnorm - s2dotnorm) * (s1dot - pdot) / (s1dot - s2dot));
        }
    }
}
