/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h>
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>

#include <time.h>

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of
	//   x, y, theta and their uncertainties from GPS) and all weights to 1.
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).

	DEBUG_PRINT1("init start.\n");
	if(0 == num_particles){
		printf("please specify num of particles.\n");
		exit(0);
	}

	default_random_engine gen;

	double std_x = std[0];
	double std_y = std[1];
	double std_theta = std[2];
	normal_distribution<double> dist_x(x, std_x);
	normal_distribution<double> dist_y(y, std_y);
	normal_distribution<double> dist_theta(theta, std_theta);

	for(int i=0; i < num_particles; i++){
		Particle p;
		p.id = i;
		p.x = dist_x(gen);
		p.y = dist_y(gen);
		p.theta = dist_theta(gen);
		p.weight = 1;
		particles.push_back(p);
		weights.push_back(p.weight);
		DEBUG_PRINT3("particle %d %f %f %f\n", particles[i].id, particles[i].x, particles[i].y, particles[i].theta);
	}
	is_initialized = true;
	DEBUG_PRINT1("init end.\n");
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

	default_random_engine gen;
	double std_x = std_pos[0];
	double std_y = std_pos[1];
	double std_theta = std_pos[2];

	DEBUG_PRINT1("prediction start.\n");

	for(int i=0; i < num_particles; i++){
		double x = particles[i].x;
		double y = particles[i].y;
		double theta = particles[i].theta;
		// Add measurements
		// CTRV Model
		if (fabs(yaw_rate) > 0.001) {
			x = x + velocity/yaw_rate * ( sin(theta + yaw_rate*delta_t) - sin(theta));
			y = y + velocity/yaw_rate * ( cos(theta) - cos(theta+yaw_rate*delta_t) );
		}
		else {
			x = x + velocity*delta_t*cos(theta);
			y = y + velocity*delta_t*sin(theta);
		}
		theta = theta + yaw_rate*delta_t;

		// Add noise (since the input is noisless, add noise here to simulate the noisy scenario???)
		normal_distribution<double> dist_x(x, std_x);
		normal_distribution<double> dist_y(y, std_y);
		normal_distribution<double> dist_theta(theta, std_theta);
		x = dist_x(gen);
		y = dist_y(gen);
		theta = dist_theta(gen);

		particles[i].x = x;
		particles[i].y = y;
		particles[i].theta = theta;
		DEBUG_PRINT3("prediction particle %d %f %f %f\n", particles[i].id, particles[i].x, particles[i].y, particles[i].theta);
	}

	DEBUG_PRINT1("prediction end.\n");
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to
	//   implement this method and use it as a helper during the updateWeights phase.

	for(int oi=0; oi < observations.size(); oi++){
		// find nearest lanmark for the observation

		if(predicted.size() == 0){
			// how to handle this situation?
			std::cout << "no lanmark around the particle!" << std::endl;
			exit(-1);
		}

		double min_dist = dist(predicted[0].x, predicted[0].y, observations[oi].x, observations[oi].y);
		observations[oi].id = predicted[0].id;

		for(int li=0; li < predicted.size(); li++){

			double dist_ = dist(predicted[li].x, predicted[li].y, observations[oi].x, observations[oi].y);

			if(dist_ <= min_dist){
				min_dist = dist_;
				observations[oi].id = predicted[li].id;
			}
		}
	}
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[],
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html

	DEBUG_PRINT1("updateWeights start.\n");

	for(int pi=0; pi < num_particles; pi++){

		Particle &particle = particles[pi];
		DEBUG_PRINT2("updateWeights particle: %d %f %f.\n", particle.id, particle.x, particle.y);

		// transform observation from vehicle coordinate to Map coordinate
		std::vector<LandmarkObs> transformed_obss;
		for(int oi=0; oi < observations.size(); oi++){
			LandmarkObs obs = observations[oi];
			// apply homogenous transformation
			LandmarkObs trans_obs;
			trans_obs.x = particle.x + cos(particle.theta)*obs.x - sin(particle.theta)*obs.y;
			trans_obs.y = particle.y + sin(particle.theta)*obs.x + cos(particle.theta)*obs.y;
			transformed_obss.push_back(trans_obs);
		}

		// Find lanmarks around current particle.
		// Does the provided observations include observations of lanmarks on the side/rear of the vehicle?
		// I don't know, so just find all lanmarks around the particle for lanmarks/observations association.
		std::vector<LandmarkObs> around_lms;
		for(int li=0; li < map_landmarks.landmark_list.size(); li ++){

			LandmarkObs lm;

			if(dist(particle.x, particle.y, \
				map_landmarks.landmark_list[li].x_f, map_landmarks.landmark_list[li].y_f) <= sensor_range)
			{
				//lm.id = map_landmarks.landmark_list[li].id_i;  // map_landmarks.landmark_list[li].id_i is start from 1 not zero
				lm.id = li;
				lm.x = map_landmarks.landmark_list[li].x_f;
				lm.y = map_landmarks.landmark_list[li].y_f;
				around_lms.push_back(lm);

				DEBUG_PRINT2("updateWeights around_lms: %d %f %f.\n", lm.id, lm.x, lm.y);
			}

		}

		dataAssociation(around_lms, transformed_obss);

		// update weights
		double weight = 1;
		for(int oi=0; oi < transformed_obss.size(); oi++){

			// Multivariate-Gaussian Probability
			double sig_x = std_landmark[0];
			double sig_y = std_landmark[1];
			double x_obs = transformed_obss[oi].x;
			double y_obs = transformed_obss[oi].y;
			double mu_x = map_landmarks.landmark_list[ transformed_obss[oi].id ].x_f;
			double mu_y = map_landmarks.landmark_list[ transformed_obss[oi].id ].y_f;

			double gauss_norm= (1./(2 * M_PI * sig_x * sig_y));
			double exponent= pow((x_obs - mu_x), 2)/(2 * pow(sig_x, 2)) + pow((y_obs - mu_y), 2)/(2 * pow(sig_y, 2));
			weight *= (gauss_norm * exp(-exponent));
			DEBUG_PRINT2("updateWeights obs:    %f %f.\n",transformed_obss[oi].x, transformed_obss[oi].y);
			DEBUG_PRINT2("updateWeights lm : %d %f %f.\n", transformed_obss[oi].id, mu_x, mu_y);
			DEBUG_PRINT2("updateWeights weight: %f .\n", weight);
		}

		particle.weight = weight;
		weights[pi] = particle.weight;
		DEBUG_PRINT2("updateWeights particle %d %f\n", particles[pi].id, particles[pi].weight);
	}
	DEBUG_PRINT1("updateWeights end.\n");
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight.
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

	DEBUG_PRINT1("resample start.\n");

	double RANDOM_MAX = 9999.0;
	std::default_random_engine gen(time(NULL));
  	std::uniform_int_distribution<int> udist(0,num_particles - 1);
	std::uniform_int_distribution<int> udist2(0,RANDOM_MAX);

	int index = udist(gen);
	double beta = 0.0;
	double mw = 0.0;

	std::vector<double>::iterator biggest = std::max_element(std::begin(weights), std::end(weights));
	mw = *biggest;

	std::vector<double> new_weights;
	std::vector<Particle> new_particles;

	// resample wheels
	for(int pi=0; pi < num_particles; pi++){

		Particle p;
		beta += 2.0 * mw * ((double)udist2(gen)/RANDOM_MAX);

		while(beta > weights[index]){
			beta -= weights[index];
			index = (index+1) % num_particles;
		}

		p = particles[index];
		new_particles.push_back(p);
		new_weights.push_back(p.weight);
	}
	weights = new_weights;
	particles = new_particles;

	DEBUG_PRINT1("resample end.\n");
}

Particle ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations,
                                     const std::vector<double>& sense_x, const std::vector<double>& sense_y)
{
    //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates

    particle.associations= associations;
    particle.sense_x = sense_x;
    particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
