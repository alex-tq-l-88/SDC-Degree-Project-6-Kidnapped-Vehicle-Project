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
#include <sstream>

#include "helper_functions.h"

using std::string;
using std::vector;

using namespace std;

#define EPS 0.00001

void ParticleFilter::init(double x, double y, double theta, double std[]) {
  /**
   * TODO: Set the number of particles. Initialize all particles to 
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1. 
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method 
   *   (and others in this file).
   */
  //Handle case that filter is already initialized
  if (is_initialized) {
    return;
  }
  
  //Initialize # of particles
  num_particles = 100;
  
  //Define normal distributions for sensor noise
  normal_distribution<double> distribution_x(x, std[0]);
  normal_distribution<double> distribution_y(y, std[1]);
  normal_distribution<double> distribution_theta(theta, std[2]);
  
  //Iniitilize particles with normal distribution with mean around value provided by GPS
  for (int i = 0; i < num_particles; i++) {
    Particle p; //create particle
    p.id = i;
    p.x = distribution_x(gen); //x dimension variation
    p.y = distribution_y(gen); //y dimension variation
    p.theta = distribution_theta(gen);
    p.weight = 1.0;

    particles.push_back(p);
  }

  // Set filter initialized status to true
  is_initialized = true;

}

void ParticleFilter::prediction(double delta_t, double std_pos[], 
                                double velocity, double yaw_rate) {
  /**
   * TODO: Add measurements to each particle and add random Gaussian noise.
   * NOTE: When adding noise you may find std::normal_distribution 
   *   and std::default_random_engine useful.
   *  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
   *  http://www.cplusplus.com/reference/random/default_random_engine/
   */
  
  //Create normal distributions
  normal_distribution<double> distribution_x(0, std_pos[0]);
  normal_distribution<double> distribution_y(0, std_pos[1]);
  normal_distribution<double> distribution_theta(0, std_pos[2]);
  
  // Update to the new state.
  for (int i = 0; i < num_particles; i++) {

  	double theta = particles[i].theta;
	
    //Check yaw rate
    if ( fabs(yaw_rate) < EPS ) { // When yaw is not changing - continues to be the same.
      particles[i].x += velocity * delta_t * cos( theta );
      particles[i].y += velocity * delta_t * sin( theta );
    } else { //update yaw rate
      particles[i].x += velocity / yaw_rate * ( sin( theta + yaw_rate * delta_t ) - sin( theta ) );
      particles[i].y += velocity / yaw_rate * ( cos( theta ) - cos( theta + yaw_rate * delta_t ) );
      particles[i].theta += yaw_rate * delta_t;
    }

    // Add noise to particle prediction
    particles[i].x += distribution_x(gen);
    particles[i].y += distribution_y(gen);
    particles[i].theta += distribution_theta(gen);
  }

}

void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted, vector<LandmarkObs>& observations) {
  /**
   * TODO: Find the predicted measurement that is closest to each 
   *   observed measurement and assign the observed measurement to this 
   *   particular landmark.
   * NOTE: this method will NOT be called by the grading code. But you will 
   *   probably find it useful to implement this method and use it as a helper 
   *   during the updateWeights phase.
   */
  
  //Iterate through each observation to find closest predicted measurement
  for (unsigned int i = 0; i < observations.size(); i++) {

    // Initialize min distance as largest possible number to ensure update
    double minDistance = numeric_limits<double>::max();

    // Initialize ID landmark - used -1 to guarantee update first time
    int mapId = -1;
    
    //Iterate through all prediction
    for (unsigned j = 0; j < predicted.size(); j++ ) { 
      //Get distance between obs & prediction in x & y dimension
      double x_Distance = observations[i].x - predicted[j].x;
      double y_Distance = observations[i].y - predicted[j].y;
      //Calculate hypothenuse
      double distance = x_Distance * x_Distance + y_Distance * y_Distance;

      //Update ID and min if pair's distance is less than prior min distance pair 
      if ( distance < minDistance ) {
        minDistance = distance;
        mapId = predicted[j].id;
      }
    }
    // Update observation's nearest prediction's ID 
    observations[i].id = mapId;
  }
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], const vector<LandmarkObs> &observations, 
                                   const Map &map_landmarks) {
  /**
   * TODO: Update the weights of each particle using a mult-variate Gaussian distribution. 
   * NOTE: The observations are given in the VEHICLE'S coordinate system. 
   *   Your particles are located according to the MAP'S coordinate system. 
   *   You will need to transform between the two systems. Keep in mind that
   *   this transformation requires both rotation AND translation (but no scaling).
   *   The following is a good resource for the theory:
   *   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
   *   and the following is a good resource for the actual equation to implement
   *   (look at equation 3.33) http://planning.cs.uiuc.edu/node99.html
   */
  
  //Constants for weight update calculations
  double std_LandmarkRange = std_landmark[0];
  double std_LandmarkBearing = std_landmark[1];

  //Iterate through each particle to updatee weights
  for (int i = 0; i < num_particles; i++) {
	//Particle's x and y dimension and direction
    double x = particles[i].x;
    double y = particles[i].y;
    double theta = particles[i].theta;
    // Find landmarks within particle's range.
    double sensor_range_2 = sensor_range * sensor_range; //sensor range squared
    vector<LandmarkObs> inRangeLandmarks; //vector of landmark observations to iterate through
    //Iterate through all landmarks
    for(unsigned int j = 0; j < map_landmarks.landmark_list.size(); j++) {
      float landmarkX = map_landmarks.landmark_list[j].x_f;
      float landmarkY = map_landmarks.landmark_list[j].y_f;
      int id = map_landmarks.landmark_list[j].id_i;
      double dX = x - landmarkX;
      double dY = y - landmarkY;
      if ( dX*dX + dY*dY <= sensor_range_2 ) { //Add landmark to in-range landmark if it is within range
        inRangeLandmarks.push_back(LandmarkObs{ id, landmarkX, landmarkY });
      }
    }

    // Transform coordinates of observation from vehicle to map dimensions
    vector<LandmarkObs> mappedObservations;
    for(unsigned int j = 0; j < observations.size(); j++) {
      double t_x = cos(theta)*observations[j].x - sin(theta)*observations[j].y + x;
      double t_y = sin(theta)*observations[j].x + cos(theta)*observations[j].y + y;
      mappedObservations.push_back(LandmarkObs{ observations[j].id, t_x, t_y });
    }

    // Use dataAssociation function to find predicted measurement closest to each observed measurement and assign observed measurement to landmark.
    dataAssociation(inRangeLandmarks, mappedObservations);

    // Reset particle weights
    particles[i].weight = 1.0;
    // Calculate particle weights.
    for(unsigned int j = 0; j < mappedObservations.size(); j++) {
      double observation_X = mappedObservations[j].x;
      double observation_Y = mappedObservations[j].y;

      int landmark_Id = mappedObservations[j].id;

      double landmark_X, landmark_Y;
      unsigned int k = 0; 
      unsigned int n_Landmarks = inRangeLandmarks.size(); //number of landmarks
      bool found = false; //no particle found by default 
      while( !found && k < n_Landmarks ) {
        if ( inRangeLandmarks[k].id == landmark_Id) {
          found = true;
          landmark_X = inRangeLandmarks[k].x;
          landmark_Y = inRangeLandmarks[k].y;
        }
        k++; //iterate K after each loop
      }

      // Calculate dx, dy and weights
      double dX = observation_X - landmark_X;
      double dY = observation_Y - landmark_Y;

      double weight = ( 1/(2*M_PI*std_LandmarkRange*std_LandmarkBearing)) * exp( -( dX*dX/(2*std_LandmarkRange*std_LandmarkRange) + (dY*dY/(2*std_LandmarkBearing*std_LandmarkBearing)) ) );
      if (weight == 0) {
        particles[i].weight *= EPS;
      } else {
        particles[i].weight *= weight;
      }
    }
  }
}

void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */
  // Get weights and max weight.
  vector<double> weights; //vector to hold weights
  double max_Weight = numeric_limits<double>::min();
  for(int i = 0; i < num_particles; i++) {
    weights.push_back(particles[i].weight);
    if ( particles[i].weight > max_Weight ) {
      max_Weight = particles[i].weight;
    }
  }

  // Initialize distributions 
  uniform_real_distribution<double> distDouble(0.0, max_Weight);
  uniform_int_distribution<int> distInt(0, num_particles - 1);

  // Initialize index and beta value
  int index = distInt(gen);
  double beta = 0.0;

  // Wheel approach demonstrated in class
  vector<Particle> resampledParticles; //vector to hold new particles
  for(int i = 0; i < num_particles; i++) {
    beta += distDouble(gen) * 2.0;  //Generate random increment
    while( beta > weights[index]) { //If beta is larger than weight of particle
      beta -= weights[index]; //Subtract particle weight from beeta
      index = (index + 1) % num_particles; //move to next particle
    }
    resampledParticles.push_back(particles[index]); //Add this particle to resampled vector
  }

  particles = resampledParticles; //update particle vector with reesample particle

}

void ParticleFilter::SetAssociations(Particle& particle, const vector<int>& associations,const vector<double>& sense_x, const vector<double>& sense_y){
  // particle: the particle to which assign each listed association, and association's (x,y) world coordinates mapping
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