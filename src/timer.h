#pragma once
#include <iostream>
#include <chrono>
#include <vector>
#include <unordered_map>

#ifndef TIMER_H
#define TIMER_H

using namespace std;

class Timer {
public:
	Timer() {

	}
	std::vector<long> times;
	std::vector<std::string> descriptions;
	std::vector<std::chrono::high_resolution_clock::time_point> starts;
	std::vector<std::chrono::high_resolution_clock::time_point> ends;

	void create(std::string name) {
		times.push_back(0);
		descriptions.push_back(name);
		//idx_map.insert(std::pair<std::string, int>(name, times.size() - 1));
		starts.push_back(std::chrono::high_resolution_clock::now());
		ends.push_back(std::chrono::high_resolution_clock::now());
	}

	void start(std::string name) {
		//int idx = idx_map[name];
		int idx = -1;
		for (int i = 0; i < descriptions.size(); i++) {
			if (descriptions[i] == name) {
				idx = i;
				break;
			}
		}
		starts[idx] = std::chrono::high_resolution_clock::now();
	}

	void end(std::string name) {
		//int idx = idx_map[name];
		int idx = -1;
		for (int i = 0; i < descriptions.size(); i++) {
			if (descriptions[i] == name) {
				idx = i;
				break;
			}
		}
		ends[idx] = std::chrono::high_resolution_clock::now();
		times[idx] +=
			std::chrono::duration_cast<std::chrono::seconds>(ends[idx] - starts[idx]).count();
	}

	void show() {
		std::cout << "Time Statistics" << std::endl;
		std::cout << std::string(20, '=') << std::endl;
		for (int i = 0; i < times.size(); i++) {
			std::cout << "Spent " << double(times[i])
				<< " seconds on " << descriptions[i] << std::endl;
		}
	}

	long log(std::string name) {
		int idx = -1;
		for (int i = 0; i < descriptions.size(); i++) {
			if (descriptions[i] == name) {
				idx = i;
				break;
			}
		}
		return times[idx];
	}
};

#endif //TIMER_H