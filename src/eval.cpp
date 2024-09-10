#include "eval.h"

//float calculate_chamfer(const m_Model& GT_model, const m_Model& Recon_model, int num_sample) {
//	std::vector<Eigen::Vector3d> vertices;
//	std::vector<Eigen::Vector3i> faces;
//	for (auto vertex : GT_model.vertices)
//		vertices.push_back(Eigen::Vector3d(vertex.x(), vertex.y(), vertex.z()));
//	for (auto face : GT_model.faces)
//		faces.push_back(Eigen::Vector3i(face.ids[0], face.ids[1], face.ids[2]));
//	open3d::geometry::TriangleMesh GT_mesh(vertices, faces);
//
//	for (auto vertex : Recon_model.vertices)
//		vertices.push_back(Eigen::Vector3d(vertex.x(), vertex.y(), vertex.z()));
//	for (auto face : Recon_model.faces)
//		faces.push_back(Eigen::Vector3i(face.ids[0], face.ids[1], face.ids[2]));
//	open3d::geometry::TriangleMesh Recon_mesh(vertices, faces);
//
//	int num_faces = faces.size();
//	if (num_faces * 10 > num_sample)
//		num_sample = num_faces * 10;
//
//	/*open3d::geometry::PointCloud gt_sample = *GT_mesh.SamplePointsPoissonDisk(num_sample);
//	open3d::geometry::PointCloud recon_sample = *Recon_mesh.SamplePointsPoissonDisk(num_sample);*/
//
//	open3d::geometry::PointCloud gt_sample = *GT_mesh.SamplePointsUniformly(num_sample);
//	open3d::geometry::PointCloud recon_sample = *Recon_mesh.SamplePointsUniformly(num_sample);
//
//	std::vector<double> dists1 = gt_sample.ComputePointCloudDistance(recon_sample);
//	std::vector<double> dists2 = recon_sample.ComputePointCloudDistance(gt_sample);
//	double sumup = 0.;
//	for (int i = 0; i < num_sample; i++)
//		sumup += dists1[i] + dists2[i];
//	return sumup / num_sample;
//}
