#include "math_util.h"


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

int python_mod(int n, int M) {
	return ((n % M) + M) % M;
}

float norm(Vector input) {
	return std::sqrt(input.squared_length());
}

std::vector<float> normalize_vector(std::vector<float>& input) {
	float length = norm(Vector(input[0], input[1], input[2]));
	for (int i = 0; i < input.size(); i++) {
		input[i] /= length;
	}
	return input;
}

Vector normalize_vector(const Vector& input) {
	float length = norm(input);
	return input / length;
}

Vector projected_vector(Vector& input, Vector& normal) {
	Vector normal_normed = normalize_vector(normal);
	return input - input * normal_normed * normal_normed;
}

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

float cal_angle_based_weight(const Vector& this_normal, const Vector& neighbor_normal) {
	float dot_pdt = abs(this_normal * neighbor_normal);
	// NOTICE!!!!!!! Possible bugs here!!!!!  can comment clamp line to see what happens on bunny.obj!!!!!
	dot_pdt = boost::algorithm::clamp(dot_pdt, 0., 1.);
	if (1. - dot_pdt < 0)
		std::cout << "error" << std::endl;
	return 1. - dot_pdt;
}

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

float cal_quality_score(Point pos_i, Point pos_u, Point pos_w,
	Vector normal_i, Vector normal_u, Vector normal_w, bool isEuclidean) {
	float len_ui = norm(pos_u - pos_i);
	float len_wi = norm(pos_w - pos_i);
	float len_uw = norm(pos_u - pos_w);
	// Project it into local 2D plane if it is too noisy.
	if (!isEuclidean) {
		Vector edge_ui = pos_u - pos_i;
		float normal_i_length = edge_ui * normalize_vector(normal_i);
		float normal_u_length = edge_ui * normalize_vector(normal_u);
		len_ui = sqrtf((len_ui * len_ui) - (normal_i_length * normal_i_length)) +
			sqrtf((len_ui * len_ui) - (normal_u_length * normal_u_length));
		len_ui /= 2.;

		Vector edge_wi = pos_w - pos_i;
		normal_i_length = edge_wi * normalize_vector(normal_i);
		float normal_w_length = edge_wi * normalize_vector(normal_w);
		len_wi = sqrtf((len_wi * len_wi) - (normal_i_length * normal_i_length)) +
			sqrtf((len_wi * len_wi) - (normal_w_length * normal_w_length));
		len_wi /= 2.;

		Vector edge_uw = pos_u - pos_w;
		normal_w_length = edge_uw * normalize_vector(normal_w);
		normal_u_length = edge_uw * normalize_vector(normal_u);
		len_uw = sqrtf((len_uw * len_uw) - (normal_w_length * normal_w_length)) +
			sqrtf((len_uw * len_uw) - (normal_u_length * normal_u_length));
		len_uw /= 2.;
	}
	float max_value = std::acos(boost::algorithm::clamp(
		(pos_w - pos_i) * (pos_u - pos_i) / len_ui /
		len_wi, -1, 1));
	float min_value = max_value;
	float radian = std::acos(boost::algorithm::clamp(
		(pos_i - pos_u) * (pos_w - pos_u) / len_ui /
		len_uw, -1, 1));
	if (radian > max_value)
		max_value = radian;
	if (radian < min_value)
		min_value = radian;
	radian = std::acos(boost::algorithm::clamp(
		(pos_u - pos_w) * (pos_i - pos_w) / len_uw /
		len_wi, -1, 1));
	if (radian > max_value)
		max_value = radian;
	if (radian < min_value)
		min_value = radian;
	return max_value - min_value;
}

float cal_edge_score(Point pos_i, Point pos_u, Point pos_w,
	Vector normal_u, Vector normal_w, float avg_edge_length, bool isEuclidean) {
	/*float len_uw = norm(pos_u - pos_w);
	if (!isEuclidean) {
		Vector edge_uw = pos_u - pos_w;
		float normal_w_length = edge_uw * normalize_vector(normal_w);
		float normal_u_length = edge_uw * normalize_vector(normal_u);
		len_uw = std::sqrtf(len_uw * len_uw - (normal_w_length * normal_w_length)) +
			std::sqrtf(len_uw * len_uw - (normal_u_length * normal_u_length));
		len_uw /= 2.;
	}
	return len_uw / (50 * avg_edge_length) * CGAL_PI;*/
	float len_ui = norm(pos_u - pos_i);
	float max_length = len_ui;
	float len_wi = norm(pos_w - pos_i);
	if (len_wi > max_length)
		max_length = len_wi;
	float len_uw = norm(pos_u - pos_w);
	if (len_uw > max_length)
		max_length = len_uw;
	return std::abs(max_length - avg_edge_length) / avg_edge_length * CGAL_PI;
}

float cal_normal_score(Vector normal_i, Vector normal_u, Vector normal_w) {
	return ((1 - normal_i * normal_u) + (1 - normal_i * normal_w)) * CGAL_PI / 4;
}

Point project_point_to_plane(Point& v, std::vector<Point>& triangle_pos,
	Vector& project_vec) {
	// result = v+kD = v+(n(p-v)/(n*D))D
	Point i = triangle_pos[0];
	Point u = triangle_pos[1];
	Point w = triangle_pos[2];
	Vector n = CGAL::cross_product(i - u, w - u);

	Point result = v + (n*(i-v)/(n*project_vec))*project_vec;
	return result;
}

bool sameSide(Point& p1, Point& p2, Point& a, Point& b) {
	Vector cp1 = CGAL::cross_product(b - a, p1 - a);
	Vector cp2 = CGAL::cross_product(b - a, p2 - a);
	if (cp1 * cp2 >= 0)
		return true;
	return false;
}

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

Vector vector2Vector(const std::vector<float>& input) {
	return Vector(input.at(0), input.at(1), input.at(2));
}

void Vector2vector(Vector input, std::vector<float>& output) {
	output.at(0) = input.x();
	output.at(1) = input.y();
	output.at(2) = input.z();
	return;
}

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

Vector triangle_mean_normal(const Vector& normal1, const Vector& normal2, const Vector& normal3) {
	Vector normal1_norm = normalize_vector(normal1);
	Vector normal2_norm = normalize_vector(normal2);
	Vector normal3_norm = normalize_vector(normal3);
	Vector output = normal1_norm + normal2_norm + normal3_norm;
	return normalize_vector(output);
}