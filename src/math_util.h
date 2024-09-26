#pragma once
#include <vector>
#include <boost/algorithm/clamp.hpp>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <Eigen/Core>
#include <Eigen/SVD>
#include <random>

#ifndef MATH_UTIL_H
#define MATH_UTIL_H

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_3 Point;
typedef Kernel::Vector_3 Vector;

void calculate_ref_vec(const Vector&, Vector&);

float norm(Vector);

double cal_radians_3d(const Vector& branch_vec, const Vector& normal);

double cal_radians_3d(const Vector&, const Vector&, const Vector& ref_vec);

float cal_angle_based_weight(const Vector&, const Vector&);

std::vector<float> normalize_vector(std::vector<float>&);

Vector normalize_vector(const Vector&);

Vector projected_vector(Vector& input, Vector& normal);

void est_normal_SVD(std::vector<int>&,
	const std::vector<Point>&, Vector&);

Point project_point_to_plane(Point&, std::vector<Point>&, Vector&);

bool sameSide(Point&, Point&, Point&, Point&);

bool point_in_triangle(Point&, std::vector<Point>&, Vector&);

Vector vector2Vector(const std::vector<float>&);

void Vector2vector(Vector, std::vector<float>&);

void add_noise(std::string, std::vector<Point>&,
	const std::vector<Vector>&, float, float, float);

void add_normal_noise(float angle, std::vector<Vector>& normals);

Vector triangle_mean_normal(const Vector&, const Vector&, const Vector&);

float cal_proj_dist(const Vector& edge,
	Vector& this_normal, Vector& neighbor_normal);

#endif
