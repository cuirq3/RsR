#include "math_util.h"

/**
 * @brief Calculate the reference vector for rotation system
 *
 * @param normal: normal direction for the target vertex
 * @param ref_vec: [OUTPUT] output the reference vector
 *
 * @return None
 */
void calculate_ref_vec(const Vector& normal, Vector& ref_vec) {
	float eps = 1e-6;
	float second = normal[1];
	if (second == 0.)
		second += eps;
	ref_vec = Vector(0, -normal[2] / second, 1);
	if (normal[2] == 1.)
		ref_vec = Vector(0., 1., 0.);
	ref_vec /= std::sqrt(ref_vec.squared_length());
	return;
}

/**
 * @brief Calculate norm of Vector
 *
 * @param input: input vector
 *
 * @return norm of Vector
 */
float norm(Vector input) {
	return std::sqrt(input.squared_length());
}

/**
 * @brief Normalize the vector to norm 1
 *
 * @param input: vector to be normalized
 *
 * @return normalized vector
 */
std::vector<float> normalize_vector(std::vector<float>& input) {
	float length = norm(Vector(input[0], input[1], input[2]));
	for (int i = 0; i < input.size(); i++) {
		input[i] /= length;
	}
	return input;
}

/**
 * @brief Normalize the Vector to norm 1
 *
 * @param input: Vector to be normalized
 *
 * @return normalized Vector
 */
Vector normalize_vector(const Vector& input) {
	float length = norm(input);
	return input / length;
}

/**
 * @brief Project a vector to a plane
 *
 * @param input: Vector to be projected
 * @param normal: normal to the plane
 *
 * @return projected Vector
 */
Vector projected_vector(Vector& input, Vector& normal) {
	Vector normal_normed = normalize_vector(normal);
	return input - input * normal_normed * normal_normed;
}

/**
 * @brief Estimate normal
 *
 * @param neighbors: index of neighbors to be used to estimate normal
 * @param vertices: the vertices of the whole point cloud
 * @param normal: [OUTPUT] estimated normal
 *
 * @return None
 */
void est_normal_SVD(std::vector<int>& neighbors,
	const std::vector<Point>& vertices, Vector& normal) {
	Eigen::Matrix<float, 3, 11> mat;
	int colidx = 0;
	for (int id : neighbors) {
		Eigen::Vector3f this_vec(vertices[id].x(), vertices[id].y(), vertices[id].z());
		mat.col(colidx) = this_vec;
		colidx++;
	}
	Eigen::Vector3f S_A = mat.rowwise().mean();
	//std::cout << mat << std::endl;
	//std::cout << S_A << std::endl;
	mat = (mat.colwise() - S_A);
	Eigen::Matrix<float, 11, 3> mat_T = mat.transpose();
	Eigen::JacobiSVD<Eigen::Matrix<float, 11, 3>> svd(mat_T, Eigen::ComputeFullU | Eigen::ComputeFullV);
	normal = Vector(svd.matrixV().transpose().row(2).x(),
		svd.matrixV().transpose().row(2).y(), svd.matrixV().transpose().row(2).z());

	//std::cout << normal.squared_length() << std::endl;
	return;
}

/**
 * @brief Calculate cos angle weight for correcting normal orientation
 *
 * @param this_normal: normal of current vertex
 * @param neighbor_normal: normal of its neighbor vertex
 *
 * @return angle weight calculated
 */
float cal_angle_based_weight(const Vector& this_normal, const Vector& neighbor_normal) {
	float dot_pdt = abs(this_normal * neighbor_normal / norm(this_normal) / norm(neighbor_normal));
	dot_pdt = boost::algorithm::clamp(dot_pdt, 0., 1.);
	if (1. - dot_pdt < 0)
		std::cout << "error" << std::endl;
	return 1. - dot_pdt;
}

/**
 * @brief Calculate the radian in the rotation system
 *
 * @param branch_vec: vector of the out-going edge
 * @param normal: normal of the root vertex
 *
 * @return radian
 */
double cal_radians_3d(const Vector& branch_vec, const Vector& normal) {
    Vector proj_vec = branch_vec - (normal * branch_vec) /
                                   norm(normal) * normal;

    Vector ref_vec = Vector(0,0,0);
    calculate_ref_vec(normal, ref_vec);

    if (norm(proj_vec) == 0.0)
        return 0.;

    Vector proj_ref = ref_vec - (normal * ref_vec) /
                                norm(normal) * normal;
    float value = boost::algorithm::clamp(
            proj_vec * proj_ref / norm(proj_vec) /
            norm(proj_ref), -1, 1);
	double radian = std::acos(value);
    if (CGAL::cross_product(proj_vec, proj_ref) * normal > 0)
        radian = 2 * CGAL_PI - radian;
    return radian;
}

/**
 * @brief Calculate the radian given the reference vector
 *
 * @param branch_vec: vector of the out-going edge
 * @param normal: normal of the root vertex
 * @param ref_vec: the reference vector
 *
 * @return radian
 */
double cal_radians_3d(const Vector& branch_vec, const Vector& normal,const Vector& ref_vec) {
	Vector proj_vec = branch_vec - (normal * branch_vec) / 
		norm(normal) * normal;
	if (norm(proj_vec) == 0.0)
		return 0.;

	Vector proj_ref = ref_vec - (normal * ref_vec) /
		norm(normal) * normal;
	float value = boost::algorithm::clamp(
		proj_vec * proj_ref / norm(proj_vec) /
		norm(proj_ref), -1, 1);
	double radian = std::acos(value);
	if (CGAL::cross_product(proj_vec, proj_ref) * normal > 0)
		radian = 2 * CGAL_PI - radian;
	return radian;
}

/**
 * @brief Project a point to the triangle plane following a given direction
 *
 * @param v: point to be projected
 * @param triangle_pos: point coordinates of 3 vertices of the triangle
 * @param project_vec: direction the point moves along
 *
 * @return the projected point
 */
Point project_point_to_plane(Point& v, std::vector<Point>& triangle_pos,
	Vector& project_vec) {

	Point i = triangle_pos[0];
	Point u = triangle_pos[1];
	Point w = triangle_pos[2];
	Vector n = CGAL::cross_product(i - u, w - u);

	Point result = v + (n * (i - v) / (n * project_vec)) * project_vec;
	return result;
}

/**
 * @brief Test if p1 and p2 are on the same side of a line defined by a and b
 *
 * @param p1: point 1
 * @param p2: point2
 * @param a: point one to define the line
 * @param b: point two to define the line
 *
 * @return true if they are on the same side
 */
bool sameSide(Point& p1, Point& p2, Point& a, Point& b) {
	Vector cp1 = CGAL::cross_product(b - a, p1 - a);
	Vector cp2 = CGAL::cross_product(b - a, p2 - a);
	if (cp1 * cp2 >= 0)
		return true;
	return false;
}

/**
 * @brief Test a point is in a triangle (2D)
 *
 * @param p: 3D points in the space
 * @param triangle_pos: point coordinates of 3 vertices of the triangle
 * @param face_normal: the normal of the local surface, different from the plane defined by the triangle
 *
 * @return true the point is inside the triangle after projection.
 */
bool point_in_triangle(Point& p, std::vector<Point>& triangle_pos,
	Vector& face_normal) {
	Point p_prime = project_point_to_plane(p, triangle_pos, face_normal);
	Point a = triangle_pos[0];
	Point b = triangle_pos[1];
	Point c = triangle_pos[2];
	if (sameSide(p_prime, a, b, c) && sameSide(p_prime, b, a, c)
		&& sameSide(p_prime, c, a, b))
		return true;
	return false;
}

/**
 * @brief Data type transfer
 *
 */
Vector vector2Vector(const std::vector<float>& input) {
	return Vector(input.at(0), input.at(1), input.at(2));
}

/**
 * @brief Data type transfer
 *
 */
void Vector2vector(Vector input, std::vector<float>& output) {
	output.at(0) = input.x();
	output.at(1) = input.y();
	output.at(2) = input.z();
	return;
}

/**
 * @brief Add noise to the point
 *
 * @param noise_type: type of the noise added
 * @param vertices: all the vertices of the point cloud
 * @param normal: corresponding normal of vertices
 * @param sigma: standard deviation of the noise
 * @param amplitude: amplitude of the noise
 * @param avg_length: average length of the kNN graph edges
 *
 * @return None
 */
void add_noise(std::string noise_type, std::vector<Point>& vertices,
	const std::vector<Vector>& normal, float sigma, float amplitude, float avg_length) {

	const int seed = 3;
	std::mt19937 mt;
	mt.seed(seed);

	std::default_random_engine generator;
	std::normal_distribution<double> normal_distribution(0, sigma);
	std::uniform_real_distribution<double> uniform_distribution(0, CGAL_PI);
	
	if (noise_type == "random") {
		for (auto& vertex : vertices) {
			float theta = uniform_distribution(mt);
			// std::cout << theta << std::endl;
			float phi = 2 * uniform_distribution(mt);
			float this_amplitude = amplitude * normal_distribution(mt);
			Vector direction(std::sin(theta)*std::cos(phi), std::sin(theta) * std::sin(phi), std::cos(theta));
			vertex += avg_length * this_amplitude * direction;
		}
	}
	else if (noise_type == "horizontal") {
		int idx = 0;
		for (auto& vertex : vertices) {
			float phi = 2 * uniform_distribution(mt);
			float this_amplitude = amplitude * normal_distribution(mt);
			Vector Ref_vec1;
			Vector this_normal = normal[idx];
			Vector direction = this_normal;
			if (norm(this_normal) != 0.) {
				calculate_ref_vec(this_normal, Ref_vec1);
				Vector Ref_vec2 = normalize_vector(CGAL::cross_product(Ref_vec1, this_normal));
				direction = Vector(std::cos(phi) * Ref_vec1 + std::sin(phi) * Ref_vec2);
			}
			if (!isfinite(norm(direction)))
				std::cout << "error" << std::endl;
			vertex += avg_length * this_amplitude * direction;
			idx++;
		}
	}
	else if (noise_type == "vertical") {
		int idx = 0;
		for (auto& vertex : vertices) {
			float this_amplitude = amplitude * normal_distribution(mt);
			Vector this_normal = normal[idx];
			Vector direction = this_normal;
			if (norm(direction) != 0.)
				direction = normalize_vector(direction);
			vertex += avg_length * this_amplitude * direction;
			idx++;
		}
	}
	else {
		std::cout << "No such noise type available now!" << std::endl;
	}
	return;
}

/**
 * @brief Add noise to the normal
 *
 * @param angle: the maximum angle the normal can pivot compared to the original direction
 * @param normals: normals exposed to noise
 *
 * @return None
 */
void add_normal_noise(float angle, std::vector<Vector>& normals) {
	const int seed = 3;
	std::mt19937 mt;
	mt.seed(seed);

	std::default_random_engine generator;
	std::uniform_real_distribution<double> uniform_distribution(0., CGAL_PI);

	for (auto& normal : normals) {
		Vector ref1;
		calculate_ref_vec(normal, ref1);
		Vector ref2 = CGAL::cross_product(normal, ref1);
		ref2 = normalize_vector(ref2);
		float phi = 2 * uniform_distribution(mt);

		Vector plane_vec = std::cos(phi) * ref1 + std::sin(phi) * ref2;
		plane_vec = normalize_vector(plane_vec);
		normal = std::cos(angle) * normal + std::sin(angle) * plane_vec;
		normal = normalize_vector(normal);
	}
}

/**
 * @brief Calculate the local surface normal (averaged direction of normals of 3 vertices in the triangle)
 *
 * @return local surface normal
 */
Vector triangle_mean_normal(const Vector& normal1, const Vector& normal2, const Vector& normal3) {
	Vector normal1_norm = normalize_vector(normal1);
	Vector normal2_norm = normalize_vector(normal2);
	Vector normal3_norm = normalize_vector(normal3);
	Vector output = normal1_norm + normal2_norm + normal3_norm;
	return normalize_vector(output);
}

/**
 * @brief Calculate projection distance
 *
 * @param edge: the edge to be considered
 * @param this_normal: normal of one vertex
 * @param neighbor_normal: normal of another vertex
 * 
 * @return projection distance
 */
float cal_proj_dist(const Vector& edge, Vector& this_normal, Vector& neighbor_normal) {
	float Euclidean_dist = norm(edge);
	float neighbor_normal_length = edge * normalize_vector(neighbor_normal);
	float normal_length = edge * normalize_vector(this_normal);
	float projection_dist = sqrtf((Euclidean_dist * Euclidean_dist) - (normal_length * normal_length));
	projection_dist += sqrtf((Euclidean_dist * Euclidean_dist) -
		(neighbor_normal_length * neighbor_normal_length));
	projection_dist /= 2.;
	if (std::abs(normalize_vector(this_normal) * normalize_vector(neighbor_normal)) < std::cos(15. / 180. * CGAL_PI))
		projection_dist = Euclidean_dist;
	return projection_dist;
}