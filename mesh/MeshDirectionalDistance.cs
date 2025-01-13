using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;

namespace g4
{
    /// <summary>
    /// Computes the minimum distance between 2 triangles, in a given direction.
    /// </summary>
    public static class MeshDirectionalDistance
    {
        public const double ContractDisplacement = 0.0001;

        public static Triangle3d ProjectTriangleOntoTrianglePlane(Triangle3d triSource, Triangle3d triTarget,
            Vector3d directionNormalized)
        {
            var triTargetPlane = new Plane3d(triTarget.V0, triTarget.V1, triTarget.V2);
            var triSourceProjPoints = new Vector3d[3];

            for (int i = 0; i < 3; ++i)
            {
                var ray = new Ray3d(triSource[i], directionNormalized, true);

                triSourceProjPoints[i] = IntrRay3Plane3.IntersectBothDirection(ray, triTargetPlane) ?? Vector3d.MaxValue;
            }

            return new Triangle3d(triSourceProjPoints[0], triSourceProjPoints[1], triSourceProjPoints[2]);
        }

        public static List<Vector3d> ComputeIntersectionPoints(Triangle3d triA, Triangle3d triB)
        {
            var intersectedPoints = new List<Vector3d>();

            var edgesTriA = GetEdges(triA);
            var edgesTriB = GetEdges(triB);

            var intersection = new IntrSegment3Segment3(new Segment3d(), new Segment3d());

            foreach (var edgeA in edgesTriA)
            {
                intersection.Segment1 = edgeA;
                foreach (var edgeB in edgesTriB)
                {
                    intersection.Segment2 = edgeB;
                    intersection.Compute();

                    if (intersection.Result != IntersectionResult.Intersects)
                    {
                        continue;
                    }

                    switch (intersection.Type)
                    {
                        case IntersectionType.Point:
                            intersectedPoints.Add(intersection.Point0);
                            break;

                        case IntersectionType.Segment:
                            intersectedPoints.Add(intersection.Point0);
                            intersectedPoints.Add(intersection.Point1);
                            break;
                    }
                }
            }

            intersectedPoints.AddRange(GetVerticesInsideTriangle(triA, triB));

            return intersectedPoints;
        }

        private static void DisplacePointsInward(List<Vector3d> points, Triangle3d triangle)
        {
            var triCenter = triangle.GetCentroid();

            for (int i = 0; i < points.Count; ++i)
            {
                points[i] += (triCenter - points[i]).Normalized * ContractDisplacement;
            }
        }

        public static bool IsCompletelyInside(Triangle3d triA, Triangle3d triB)
        {
            return triB.IsPointInTriangle(triA.V0) &&
                   triB.IsPointInTriangle(triA.V1) &&
                   triB.IsPointInTriangle(triA.V2);
        }

        public static List<Vector3d> GetVerticesInsideTriangle(Triangle3d triA, Triangle3d triB)
        {
            var pointsInside = new List<Vector3d>(3);

            for (int i = 0; i < 3; ++i)
            {
                if (triB.IsPointInTriangle(triA[i]))
                {
                    pointsInside.Add(triA[i]);
                }
            }

            return pointsInside;
        }

        public static (double, Segment3d) GetMinimumDistanceIntersection(Ray3d[] rays, Triangle3d triangle)
        {
            var minDist = double.PositiveInfinity;
            var closetOrigin = Vector3d.MaxValue;
            var closestHitPoint = Vector3d.MaxValue;

            var rayInt = new IntrRay3Triangle3(new Ray3d(), new Triangle3d());

            for (int i = 0; i < rays.Length; ++i)
            {
                var ray = rays[i];

                rayInt.Ray = rays[i];
                rayInt.Triangle = triangle;
                rayInt.Compute();

                if (rayInt.Result == IntersectionResult.Intersects)
                {
                    var hitPoint = ray.PointAt(rayInt.RayParameter);
                    var rayHitDist = ray.Origin.Distance(hitPoint);

                    if (rayHitDist < minDist)
                    {
                        minDist = rayHitDist;
                        closetOrigin = ray.Origin;
                        closestHitPoint = hitPoint;
                    }
                }
            }

            return (minDist, new Segment3d(closetOrigin, closestHitPoint));
        }

        private static Segment3d[] GetEdges(Triangle3d triangle)
        {
            return new[]
            {
                new Segment3d(triangle.V0, triangle.V1),
                new Segment3d(triangle.V1, triangle.V2),
                new Segment3d(triangle.V2, triangle.V0)
            };
        }

        public static (double distance, Vector3d point1, Vector3d point2) ComputeMinimumDirectionalDistance(
            Triangle3d tri1, Triangle3d tri2, Vector3d directionNormalized)
        {
            Triangle3d tri1Proj = ProjectTriangleOntoTrianglePlane(tri1, tri2, directionNormalized);

            if (tri1Proj.V0 == Vector3d.MaxValue ||
                tri1Proj.V1 == Vector3d.MaxValue ||
                tri1Proj.V2 == Vector3d.MaxValue)
            {
                return (double.PositiveInfinity, Vector3d.MaxValue, Vector3d.MaxValue);
            }

            var intersectionPoints = ComputeIntersectionPoints(tri1Proj, tri2);

            // Contracts slightly inwards the triangle, to avoid points exactly on the edge
            if (IsCompletelyInside(tri1Proj, tri2))
            {
                DisplacePointsInward(intersectionPoints, tri1Proj);
            }
            else
            {
                // Points are on the boundary of tri2
                DisplacePointsInward(intersectionPoints, tri2);
            }

            // Intersection points + points from t2 inside t1 projection
            intersectionPoints.AddRange(GetVerticesInsideTriangle(tri2, tri1Proj));

            // Cast rays back and take the shortest distance
            var rays = intersectionPoints.Select(p => new Ray3d(p, -directionNormalized, true)).ToArray();
            (double closesDist, Segment3d closestSeg) = GetMinimumDistanceIntersection(rays, tri1);

            // Found minimum distance
            if (closestSeg.P0 != Vector3d.MaxValue && closestSeg.P1 != Vector3d.MaxValue)
            {
                return (closesDist, closestSeg.P0, closestSeg.P1);
            }

            return (double.PositiveInfinity, Vector3d.MaxValue, Vector3d.MaxValue);
        }

        public static (double distance, Vector3d point1, Vector3d point2) ComputeMinimumDirectionalDistance(
            DMesh3 meshA, DMesh3 meshB, Vector3d directionNormalized)
        {
            var triAIndices = meshA.TriangleIndices().ToArray();
            var triBIndices = meshB.TriangleIndices().ToArray();

            (double distance, Vector3d point1, Vector3d point2) result = (double.PositiveInfinity, Vector3d.MaxValue, Vector3d.MaxValue);

            var syncObj = new object();

            Parallel.ForEach(triAIndices, triAIndex =>
            {
                Triangle3d triA = default;
                meshA.GetTriVertices(triAIndex, ref triA.V0, ref triA.V1, ref triA.V2);

                foreach (var tri2Index in triBIndices)
                {
                    Triangle3d triB = default;
                    meshB.GetTriVertices(tri2Index, ref triB.V0, ref triB.V1, ref triB.V2);

                    if (triA.IsDegenerate() ||
                        triB.IsDegenerate())
                    {
                        continue;
                    }

                    var minDirectionalDistance = ComputeMinimumDirectionalDistance(triA, triB, directionNormalized);

                    lock (syncObj)
                    {
                        if (minDirectionalDistance.distance < result.distance)
                        {
                            result = minDirectionalDistance;
                        }
                    }
                }
            });

            return result;
        }
    }
}