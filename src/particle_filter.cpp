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

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of
	//   x, y, theta and their uncertainties from GPS) and all weights to 1.
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
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
		DEBUG_PRINT("particle %d %f %f %f\n", particles[i].id, particles[i].x, particles[i].y, particles[i].theta);
	}
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
	}
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to
	//   implement this method and use it as a helper during the updateWeights phase.

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

	for(int pi=0; pi < num_particles; pi++){

		Particle &particle = particles[pi];

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
				lm.id = map_landmarks.landmark_list[li].id_i;
				lm.x = map_landmarks.landmark_list[li].x_f;
				lm.y = map_landmarks.landmark_list[li].y_f;
				around_lms.push_back(lm);
			}
		}

		dataAssociation(around_lms, transformed_obss);




	}
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight.
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

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
