/**
 * particle_filter.cpp
 *
 * Created on: Dec 12, 2016
 * Author: Tiffany Huang
 */

#include "particle_filter.h"

#include <math.h>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>
#include <string>
#include <vector>

#include "helper_functions.h"

using std::string;
using std::vector;
using namespace std;

// declare a random engine to be used across multiple and various method calls
static default_random_engine gen;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
  /**
   * TODO: Set the number of particles. Initialize all particles to 
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1. 
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method 
   *   (and others in this file).
   */
  num_particles = 100;  // TODO: Set the number of particles
  
  normal_distribution<double> dist_x(x, std[0]);
  normal_distribution<double> dist_y(y, std[1]);
  normal_distribution<double> dist_theta(theta, std[2]);
    
  for (int i = 0; i < num_particles; i++){
    Particle particle;
    particle.id 	= i;
    particle.x 		= dist_x(gen);
    particle.y 		= dist_y(gen);
    particle.theta 	= dist_theta(gen);
    particle.weight	= 1.0;
    
    particles.push_back(particle);
    }
  
  is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], 
                                double velocity, double yaw_rate){
  /**
   * TODO: Add measurements to each particle and add random Gaussian noise.
   * NOTE: When adding noise you may find std::normal_distribution 
   *   and std::default_random_engine useful.
   *  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
   *  http://www.cplusplus.com/reference/random/default_random_engine/
   */
  normal_distribution<double> dist_x(0, std_pos[0]);
  normal_distribution<double> dist_y(0, std_pos[1]);
  normal_distribution<double> dist_theta(0, std_pos[2]);
  
// Add measurements to particles, Different equations depending on yaw_rate = zero or not
  for (int i = 0; i < num_particles; i++){
    if (yaw_rate == 0){
      particles[i].x		+= velocity * delta_t * cos(particles[i].theta);
      particles[i].y		+= velocity * delta_t * sin(particles[i].theta);
      particles[i].theta	+= 0; // Theta will stay the same due to no yaw_rate
    }
    else{
      particles[i].x 		+= (velocity/yaw_rate) * (sin(particles[i].theta + (yaw_rate * delta_t)) - sin(particles[i].theta));
      particles[i].y 		+= (velocity/yaw_rate) * (cos(particles[i].theta) - cos(particles[i].theta + (yaw_rate * delta_t)));
      particles[i].theta	+= yaw_rate * delta_t;      
    }

    // Add noise to the particles
    particles[i].x		+= dist_x(gen);
    particles[i].y		+= dist_y(gen);
    particles[i].theta	+= dist_theta(gen);
  }
}

void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted, 
                                     vector<LandmarkObs>& observations) {
  /**
   * TODO: Find the predicted measurement that is closest to each 
   *   observed measurement and assign the observed measurement to this 
   *   particular landmark.
   * NOTE: this method will NOT be called by the grading code. But you will 
   *   probably find it useful to implement this method and use it as a helper 
   *   during the updateWeights phase.
   */
 for (unsigned int i = 0; i < observations.size(); i++){    
    // init minimum distance to maximum possible, id of landmark from map placeholder
    double min_dist = numeric_limits<double>::max();
    int map_id = -1;
    
    for (unsigned int j = 0; j < predicted.size(); j++){
      // get distance between current/predicted landmarks
      double cur_dist = dist(observations[i].x, observations[i].y, predicted[j].x, predicted[j].y);

      // find the predicted landmark nearest the current observed landmark
      if (cur_dist < min_dist) {
        min_dist = cur_dist;
        map_id = predicted[j].id;
      }
    }

    // set the observation's id to the nearest predicted landmark's id
    observations[i].id = map_id;
  }
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
                                   const vector<LandmarkObs> &observations, 
                                   const Map &map_landmarks) {
  /**
   * TODO: Update the weights of each particle using a mult-variate Gaussian 
   *   distribution. You can read more about this distribution here: 
   *   https://en.wikipedia.org/wiki/Multivariate_normal_distribution
   * NOTE: The observations are given in the VEHICLE'S coordinate system. 
   *   Your particles are located according to the MAP'S coordinate system. 
   *   You will need to transform between the two systems. Keep in mind that
   *   this transformation requires both rotation AND translation (but no scaling).
   *   The following is a good resource for the theory:
   *   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
   *   and the following is a good resource for the actual equation to implement
   *   (look at equation 3.33) http://planning.cs.uiuc.edu/node99.html
   */
// Process for each particle...
  for (int i = 0; i < num_particles; i++){

    // get the particle x, y coordinates and bearing theta
    double p_x = particles[i].x;
    double p_y = particles[i].y;
    double p_theta = particles[i].theta;

    // create a vector(transformed_obs) of observations, which are transformed from vehicle coordinates to map coordinates
    vector<LandmarkObs> transformed_obs;
    
    // apply Homogeneous transformation equation to each observations
    for (unsigned int j = 0; j < observations.size(); j++) {
      double tobs_x = cos(p_theta)*observations[j].x - sin(p_theta)*observations[j].y + p_x;
      double tobs_y = sin(p_theta)*observations[j].x + cos(p_theta)*observations[j].y + p_y;
      transformed_obs.push_back(LandmarkObs{ observations[j].id, tobs_x, tobs_y });
    }

    // Create a vector(sensored_landmarks) to hold the map landmark locations within sensor range of each particle
    vector<LandmarkObs> sensored_landmarks;

    // Process for each landmark...
    for (unsigned int j = 0; j < map_landmarks.landmark_list.size(); j++) {

      // get id and x,y coordinates
      int lm_id = map_landmarks.landmark_list[j].id_i;
      float lm_x = map_landmarks.landmark_list[j].x_f;
      float lm_y = map_landmarks.landmark_list[j].y_f;
      
      // only consider landmarks within sensor range of the particle
      if ( dist(lm_x,lm_y,p_x,p_y) <= sensor_range ) {
        sensored_landmarks.push_back(LandmarkObs{ lm_id, lm_x, lm_y });
      }
    }


    // perform dataAssociation for the predictions and transformed observations on current particle
    dataAssociation(sensored_landmarks, transformed_obs);

    // initialize particle weight for each time stamp
    particles[i].weight = 1.0;

    for (unsigned int j = 0; j < transformed_obs.size(); j++) {
      
      // placeholders for observation and associated prediction coordinates
      double o_x, o_y, pr_x, pr_y;
      o_x = transformed_obs[j].x;
      o_y = transformed_obs[j].y;
      int associated_prediction = transformed_obs[j].id;

      // get the x,y coordinates of the prediction associated with the current observation
      for (unsigned int k = 0; k < sensored_landmarks.size(); k++) {
        if (sensored_landmarks[k].id == associated_prediction) {
          pr_x = sensored_landmarks[k].x;
          pr_y = sensored_landmarks[k].y;
        }
      }

      // calculate weight for this observation with multivariate Gaussian
      double s_x = std_landmark[0];
      double s_y = std_landmark[1];
      double obs_w = ( 1/(2*M_PI*s_x*s_y)) * exp( -( pow(pr_x-o_x,2)/(2*pow(s_x, 2)) + (pow(pr_y-o_y,2)/(2*pow(s_y, 2))) ) );

      // product of this obersvation weight with total observations weight
      particles[i].weight *= obs_w;
    }
    
  }  // end of Process for each particle...
}

void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */
  vector<Particle> new_particles;

  // get all of the current weights
  vector<double> weights;
  for (int i = 0; i < num_particles; i++) {
    weights.push_back(particles[i].weight);
  }

  // generate random starting index for resampling wheel
  uniform_int_distribution<int> uniintdist(0, num_particles-1);
  auto index = uniintdist(gen);

  // get max weight
  double max_weight = *max_element(weights.begin(), weights.end());

  // uniform random distribution [0.0, max_weight)
  uniform_real_distribution<double> unirealdist(0.0, max_weight);

  double beta = 0.0;

  // spin the resample wheel!
  for (int i = 0; i < num_particles; i++) {
    beta += unirealdist(gen) * 2.0;
    while (beta > weights[index]) {
      beta -= weights[index];
      index = (index + 1) % num_particles;
    }
    new_particles.push_back(particles[index]);
  }

  particles = new_particles;
}

void ParticleFilter::SetAssociations(Particle& particle, 
                                     const vector<int>& associations, 
                                     const vector<double>& sense_x, 
                                     const vector<double>& sense_y) {
  // particle: the particle to which assign each listed association, 
  //   and association's (x,y) world coordinates mapping
  // associations: The landmark id that goes along with each listed association
  // sense_x: the associations x mapping already converted to world coordinates
  // sense_y: the associations y mapping already converted to world coordinates
  particle.associations= associations;
  particle.sense_x = sense_x;
  particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best) {
  vector<int> v = best.associations;
  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<int>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}

string ParticleFilter::getSenseCoord(Particle best, string coord) {
  vector<double> v;

  if (coord == "X") {
    v = best.sense_x;
  } else {
    v = best.sense_y;
  }

  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}