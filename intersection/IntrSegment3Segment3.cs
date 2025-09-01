using System;

namespace g4
{
    public class IntrSegment3Segment3
    {
        private Segment3d _segment1;

        private Segment3d _segment2;

        private double _epsilon = 1E-08;

        public IntersectionResult Result;

        public IntersectionType Type;

        public Vector3d Point0;

        public Vector3d Point1;

        public Segment3d Segment1
        {
            get => _segment1;
            set
            {
                _segment1 = value;
                Result = IntersectionResult.NotComputed;
            }
        }

        public Segment3d Segment2
        {
            get => _segment2;
            set
            {
                _segment2 = value;
                Result = IntersectionResult.NotComputed;
            }
        }

        public double Epsilon
        {
            get => _epsilon;
            set
            {
                _epsilon = Math.Max(value, 0.0);
                Result = IntersectionResult.NotComputed;
            }
        }

        public IntrSegment3Segment3(Segment3d seg1, Segment3d seg2)
        {
            _segment1 = seg1;
            _segment2 = seg2;

            Result = IntersectionResult.NotComputed;
            Type = IntersectionType.Empty;
        }

        public IntrSegment3Segment3 Compute()
        {
            Find();
            return this;
        }

        private bool HandleCollinearSegments(Segment3d seg1, Segment3d seg2)
        {
            var p1 = seg1.P0;
            var q1 = seg1.P1;
            var p2 = seg2.P0;
            var q2 = seg2.P1;

            // Sort points by projection along segment direction
            var t0 = Project(p1, q1, p2);
            var t1 = Project(p1, q1, q2);

            if (t0 > t1) {
	            var tt = t0;
	            t0 = t1;
	            t1 = tt;
            }

            // Check overlap
            if (t1 < 0.0 || t0 > 1.0)
            {
                Result = IntersectionResult.NoIntersection;
                Type = IntersectionType.Empty;

                return false;
            }

            // Return the overlapping segment
            var overlapStart = p1 + Math.Max(0, t0) * (q1 - p1);
            var overlapEnd = p1 + Math.Min(1, t1) * (q1 - p1);

            if (overlapStart.Distance(overlapEnd) < _epsilon)
            {
                Result = IntersectionResult.Intersects;
                Type = IntersectionType.Point;
                Point0 = overlapStart;

                return true;
            }

            Result = IntersectionResult.Intersects;
            Type = IntersectionType.Segment;
            Point0 = overlapStart;
            Point1 = overlapEnd;

            return true;
        }

        private double Project(Vector3d p1, Vector3d q1, Vector3d point)
        {
            var direction = q1 - p1;
            var squaredLength = direction.LengthSquared;

            if (squaredLength < _epsilon)
            {
                return 0.0;
            }

            return (point - p1).Dot(direction) / squaredLength;
        }

        public bool Find()
        {
            var p1 = _segment1.P0;
            var q1 = _segment1.P1;
            var p2 = _segment2.P0;
            var q2 = _segment2.P1;

            var d1 = q1 - p1; // First segment direction
            var d2 = q2 - p2; // Second segment direction

            // Solve system to find t and u
            var r = p2 - p1;
            var crossD1D2 = d1.Cross(d2);
            var denominator = crossD1D2.LengthSquared;

            if (denominator < _epsilon)
            {
                // Check if they are collinear
                if (r.Cross(d1).LengthSquared < _epsilon)
                {
                    // Collinear: check overlap
                    return HandleCollinearSegments(_segment1, _segment2);
                }

                Result = IntersectionResult.NoIntersection;
                Type = IntersectionType.Empty;
                return false; // Non-collinear parallels
            }

            // Compute t and u parameters
            var t = r.Cross(d2).Dot(crossD1D2) / denominator;
            var u = r.Cross(d1).Dot(crossD1D2) / denominator;

            // Check if t and u are in the range [0, 1]
            if (t < 0.0 || t > 1.0 || u < 0.0 || u > 1.0)
            {
                // Intersection outside segments
                Result = IntersectionResult.NoIntersection;
                Type = IntersectionType.Empty;
                return false;
            }

            var intersectionPoint = p1 + t * d1;
            Result = IntersectionResult.Intersects;
            Type = IntersectionType.Point;
            Point0 = intersectionPoint;

            return true;
        }
    }
}

