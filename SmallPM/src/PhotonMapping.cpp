/*********************************************************************************
Copyright (C) 2014 Adrian Jarabo (ajarabo@unizar.es)
Copyright (C) 2014 Diego Gutierrez (diegog@unizar.es)
All rights reserved.

This is an educational Ray Tracer developed for the course 'Informatica Grafica'
(Computer Graphics) tought at Universidad de Zaragoza (Spain). As such, it does not 
intend to be fast or general, but just to provide an educational tool for undergraduate
students. 

This software is provided as is, and any express or implied warranties are disclaimed.
In no event shall copyright holders be liable for any damage.
**********************************************************************************/
#include "PhotonMapping.h"
#include "World.h"
#include "Intersection.h"
#include "Ray.h"
#include "BSDF.h"

const int MAX_PHOTONS = 100;

//*********************************************************************
// Compute the photons by tracing the Ray 'r' from the light source
// through the scene, and by storing the intersections with matter
// in the lists 'xx_photons', storing the diffuse (global) and caustic
// photons respectively. For efficiency, both are computed at the same
// time, since computing them separately would result into a lost of
// several samples marked as caustic or diffuse.
// Same goes with the boolean 'direct', that specifies if direct 
// photons (from light to surface) are being stored or not. 
// The initial traced photon has energy defined by the tristimulus
// 'p', that accounts for the emitted power of the light source.
// The function will return true when there are more photons (caustic
// or diffuse) to be shot, and false otherwise.
//---------------------------------------------------------------------
bool PhotonMapping::trace_ray(const Ray& r, const Vector3 &p, 
			   std::list<Photon> &global_photons, std::list<Photon> &caustic_photons, bool direct)
{

	//Check if max number of shots done...
	if( ++m_nb_current_shots > m_max_nb_shots )
	{
		return false;
	}
	
	// Compute irradiance photon's energy
	Vector3 energy(p);
	
	Ray photon_ray(r);
	photon_ray.shift();

	bool is_caustic_particle = false;

	//Iterate the path
	while(1)
	{
		// Throw ray and update current_it
		Intersection it;
		world->first_intersection(photon_ray, it);

		if( !it.did_hit() )
			break;

		//Check if has hit a delta material...
		if( it.intersected()->material()->is_delta() )
		{
			// If delta material, then is caustic...
			// Don't store the photon!
			is_caustic_particle = true;
		}
		else if (photon_ray.get_level() > 0 || direct)
		{
			//If non-delta material, store the photon!
			if( is_caustic_particle )	
			{				
				//If caustic particle, store in caustics
				if( caustic_photons.size() < m_nb_caustic_photons )
					caustic_photons.push_back( Photon(it.get_position(), photon_ray.get_direction(), energy ));
			}
			else						
			{
				//If non-caustic particle, store in global
				if( global_photons.size() < m_nb_global_photons )
					global_photons.push_back( Photon(it.get_position(), photon_ray.get_direction(), energy ));
			}
			is_caustic_particle = false;
		}	
		
		Real pdf;

		Vector3 surf_albedo = it.intersected()->material()->get_albedo(it);
		Real avg_surf_albedo = surf_albedo.avg();

		Real epsilon2 = static_cast<Real>(rand())/static_cast<Real>(RAND_MAX);
		while (epsilon2 < 0.)
			epsilon2 = static_cast<Real>(rand())/static_cast<Real>(RAND_MAX);
		
		if (epsilon2 > avg_surf_albedo || photon_ray.get_level() > 20 ) 
			break;
			
		// Random walk's next step
		// Get sampled direction plus pdf, and update attenuation
		it.intersected()->material()->get_outgoing_sample_ray(it, photon_ray, pdf );

		// Shade...
		energy = energy*surf_albedo;
		if( !it.intersected()->material()->is_delta() )
			energy *= dot_abs(it.get_normal(), photon_ray.get_direction())/3.14159;		

		energy = energy /(pdf*avg_surf_albedo);
	}
	
	if( caustic_photons.size() == m_nb_caustic_photons && 
		global_photons.size() == m_nb_global_photons )
	{
		m_max_nb_shots = m_nb_current_shots-1;
		return false;
	}

	return true;
}

//*********************************************************************
// TODO: Implement the preprocess step of photon mapping,
// where the photons are traced through the scene. To do it,
// you need to follow these steps for each shoot:
//  1 - Sample a world's light source in the scene to create
//		the initial direct photon from the light source.
//	2 - Trace the photon through the scene storing the inter-
//		sections between the photons and matter. You can use
//		the function 'trace_ray' for this purpose.
//	3 - Finally, once all the photons have been shot, you'll
//		need to build the photon maps, that will be used later
//		for rendering. 
//		NOTE: Careful with function
//---------------------------------------------------------------------
void PhotonMapping::preprocess()
{
	
	std::vector<LightSource*> ls = world->light_source_list;
	Vector3 centro(0, 0, 0);
	int radio = 2;
	double fmin = -(radio / 2);
	std::list<Photon> *global_photons = new std::list<Photon>();
	std::list<Photon> *caustic_photons = new std::list<Photon>();
	Vector3 flujo(1, 1, 1);
	Ray *r;
	for (size_t i = 0; i < ls.size(); i++)
	{
		LightSource *luz = ls[i];
		do{
			// Hace lo del cuadradico
			// Create random point with rante from -1 to 1

			double x ;
			double y ;
			double z ;
			double dist = radio+1;
			// Itera hasta que encuentra un punto que este dentro de la esfera
			while (!(dist<= radio)){
				x = fmin + ((double)rand() / RAND_MAX) * radio;
				y = fmin + ((double)rand() / RAND_MAX) * radio;
				z = fmin + ((double)rand() / RAND_MAX) * radio;				
				dist = sqrt(x*x + y*y + z*z);
			}
			Vector3* random = new Vector3(x, y, z);
			r = new Ray(luz->get_position(), centro - *random, 0);

		} while (trace_ray(*r, flujo, *global_photons, *caustic_photons, false));
		unsigned int n_phot_g = m_max_nb_shots;
		int n_global = global_photons->size();
		for (int i = 0; i < n_global; i++){
			Photon photon = global_photons->front();
			global_photons->pop_front();
			m_global_map.store(std::vector<Real>(photon.position.data,
				photon.position.data + 3), photon);
			m_global_map.balance();
		}

		int n_caustics = caustic_photons->size();
		for (int i = 0; i < n_caustics; i++){
			Photon photon = caustic_photons->front();
			caustic_photons->pop_front();
			m_caustics_map.store(std::vector<Real>(photon.position.data,
				photon.position.data + 3), photon);
			m_caustics_map.balance();
		}
	}
	
}


double PhotonMapping::distancia(Vector3* a, Vector3 b){
	return sqrt(pow(a->getComponent(1) - b.getComponent(1), 2) +
		pow(a->getComponent(2) - b.getComponent(2), 2) +
		pow(a->getComponent(3) - b.getComponent(3), 2)
		);
}

//*********************************************************************
// TODO: Implement the function that computes the rendering equation 
// using radiance estimation with photon mapping, using the photon
// maps computed as a proprocess. Note that you will need to handle
// both direct and global illumination, together with recursive the 
// recursive evaluation of delta materials. For an optimal implemen-
// tation you should be able to do it iteratively.
// In principle, the class is prepared to perform radiance estimation
// using k-nearest neighbors ('m_nb_photons') to define the bandwidth
// of the kernel.
//---------------------------------------------------------------------
Vector3 PhotonMapping::shade(Intersection &it0)const
{
	Vector3 L(0);
	Intersection it(it0);

	// Calculates direct light
	Vector3 kd = it.intersected()->material()->get_albedo(it);
	Vector3 normal = normalize(it.get_normal());
	std::vector<LightSource*> world_lights = world->light_source_list;
	LightSource* ls = world_lights[0]; // Falta tener en cuenta todas las luces
	Vector3 light = normalize(ls->get_position() - it.get_position());
	double escProd = dot(normal, light);
	// If angle < 0, it is not illuminated so intensity = 0
	if (escProd < 0) {
		escProd = 0;
	}
	Vector3 intensity = kd * escProd;

	// Photon mapping algorithm for Global Illumination
	std::vector<const KDTree<Photon, 3>::Node*> global_photons;
	Real max_distance = 1;
	do {
		m_global_map.find(std::vector<Real>(it.get_position, it.get_position + 3),
			m_nb_photons, global_photons, max_distance);
		max_distance++;
	} while (global_photons.size < MAX_PHOTONS);

	// Photon mapping algorithm for Caustics
	std::vector<const KDTree<Photon, 3>::Node*> caustic_photons;
	max_distance = 1;
	do {
		m_global_map.find(std::vector<Real>(it.get_position, it.get_position + 3),
			m_nb_photons, caustic_photons, max_distance);
		max_distance++;
	} while (caustic_photons.size < MAX_PHOTONS);

	return intensity;

	//**********************************************************************
	// The following piece of code is included here for two reasons: first
	// it works as a 'hello world' code to check that everthing compiles 
	// just fine, and second, to illustrate some of the functions that you 
	// will need when doing the work. Goes without saying: remove the 
	// pieces of code that you won't be using.
	//
	unsigned int debug_mode = 1;

	switch (debug_mode)
	{
	case 1:
		// ----------------------------------------------------------------
		// Display Albedo Only
		L = it.intersected()->material()->get_albedo(it);
		break;
	case 2:
		// ----------------------------------------------------------------
		// Display Normal Buffer
		L = it.get_normal();
		break;
	case 3:
		// ----------------------------------------------------------------
		// Display whether the material is specular (or refractive) 
		L = Vector3(it.intersected()->material()->is_delta());
		break;

	case 4:
		// ----------------------------------------------------------------
		// Display incoming illumination from light(0)
		L = world->light(0).get_incoming_light(it.get_position());
		break;

	case 5:
		// ----------------------------------------------------------------
		// Display incoming direction from light(0)
		L = world->light(0).get_incoming_direction(it.get_position());
		break;

	case 6:
		// ----------------------------------------------------------------
		// Check Visibility from light(0)
		if (world->light(0).is_visible(it.get_position()))
			L = Vector3(1.);
		break;
	}
	// End of exampled code
	//**********************************************************************

	return L;
}