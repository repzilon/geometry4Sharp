﻿namespace g4
{
    public class DCurveProjectionTarget : IProjectionTarget
    {
        public DCurve3 Curve;


        public DCurveProjectionTarget(DCurve3 curve)
        {
            Curve = curve;
        }

        public Vector3d Project(Vector3d vPoint, int identifier = -1)
        {
            Vector3d vNearest = Vector3d.Zero;
            double fNearestSqr = double.MaxValue;

            int N = Curve.VertexCount;
            int NStop = (Curve.Closed) ? N : N - 1;
            for ( int i = 0; i < NStop; ++i ) {
                Segment3d seg = new Segment3d(Curve[i], Curve[(i + 1) % N]);
                Vector3d pt = seg.NearestPoint(vPoint);
                double dsqr = pt.DistanceSquared(vPoint);
                if (dsqr < fNearestSqr) {
                    fNearestSqr = dsqr;
                    vNearest = pt;
                }
            }

            return (fNearestSqr < double.MaxValue) ? vNearest : vPoint;
        }
    }
}
