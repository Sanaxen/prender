#ifndef _CIRCLE_H
#define _CIRCLE_H

#include "entity.h"
#include "plane.h"

#include"KerrBlackHole.h"
#include "WormHole.h"

namespace prender {
	extern KerrBlackHole* BlackHole;

	class DiscMapping
	{
	private:
		double rMin;
		double rMax;
		int sizeX;
		int sizeY;

	public:

		DiscMapping(double rMin_, double rMax_, int sizex_, int sizey_)
		{
			rMax = rMax_;
			rMin = rMin_;
			sizeX = sizex_;
			sizeY = sizey_;
		}


		void Map(double r, double theta, double phi, int& x, int& y)
		{
			//y = sizeY-1 - y;
			if (r < rMin || r > rMax)
			{
				x = 0;
				y = sizeY;
			}

			x = (int)(phi / (2 * PS_PI) * sizeX) % sizeX;
			if (x < 0) x = sizeX + x;
			y = (int)((r - rMin) / (rMax - rMin) * sizeY);
			if (y > sizeY - 1)
				y = sizeY - 1;
			if (y < 0) y = 0;
			//y = sizeY - y - 1;

		}
		void MapUV(double disk_r, double r, double u, double v, int& x, int& y)
		{
			//y = sizeY-1 - y;
			if (r < rMin || r > disk_r)
			{
				x = -1;
				y = -1;
			}

			double theta = atan2(u/r, v/r);

			r = (r - rMin) / (disk_r - rMin);
			if (r < 0)
			{
				x = -1;
				y = -1;
			}
			x = (sizeX*0.5 + sizeX*r*cos(theta));
			y = (sizeY*0.5 + sizeY*r*sin(theta));
			if (x < 0 || x >= sizeX || y < 0 || y >= sizeY)
			{
				x = -1;
				y = -1;
			}
		}
	};

	class SphericalMapping
	{
		int SizeX;
		int SizeY;

	public:
		SphericalMapping(int sizex_, int sizey_)
		{
			SizeX = sizex_;
			SizeY = sizey_;
		}

		void Map(double r, double theta, double phi, int& x, int& y)
		{
			// do mapping of texture image
			double textureScale = 1.0;

			x = (int)(((phi * textureScale) / (2 * PS_PI)) * this->SizeX) % this->SizeX;
			y = (int)((theta * textureScale / PS_PI) * this->SizeY) % this->SizeY;

			if (x < 0) x = this->SizeX + x;
			if (y < 0) y = this->SizeY + y;

			//y = SizeY - y - 1;

		}
	};

	extern std::vector<DiscMapping*> discMap;
	extern std::vector<BitMap*> discTexture;

	class Circle : public UVPlane
	{
	public:
		UVPlane* base;
		double radius;
		double radius_inner = 0.0;
		int blackhole_disk;

		Circle(const Vector3d& o, const Vector3d& n, const double r, int m)
		{
			type = EntityType::ENTITY_TYPE_CIRCLE;
			material_id = m;
			org = o;

			blackhole_disk = -1;
			Vector3d orienting_normal, w, u, v;
			w = n;
			orienting_normal = n;
			OrthonormalBasis(orienting_normal, w, u, v);

			radius = r;
			base = new UVPlane(o, r*u, r*v, m);
			base->circle = true;
		}

		~Circle()
		{
			delete base;
		}
		void CalcArea()
		{
			area = radius*radius*PS_PI;
		}

		Entity* ConstructBVH()
		{
			return this;
		}
		Entity* ConstructQBVH()
		{
			return this;
		}
		void CreateBoundingBox()
		{
			const double s1 = -1.0;
			const double s2 = 1.0;
			Vector3d q[4];
			q[0] = org + u_length*s1*u_axis + v_length*s1*v_axis;
			q[1] = org + u_length*s1*u_axis + v_length*s2*v_axis;
			q[2] = org + u_length*s2*u_axis + v_length*s1*v_axis;
			q[3] = org + u_length*s2*u_axis + v_length*s2*v_axis;

			Vector3d min = q[0];
			Vector3d max = q[0];
			for (int i = 1; i < 4; i++)
			{
				if (min.x > q[i].x) min.x = q[i].x;
				if (min.y > q[i].y) min.y = q[i].y;
				if (min.z > q[i].z) min.z = q[i].z;
				if (max.x < q[i].x) max.x = q[i].x;
				if (max.y < q[i].y) max.y = q[i].y;
				if (max.z < q[i].z) max.z = q[i].z;
			}
			boundingBox = BoundingBox(min, max);
		}

		bool intersect(const Ray &ray, IntersectionPos *hitpoint) const
		{
			if (!base->intersect(ray, hitpoint))
			{
				return false;
			}
			if (dot(hitpoint->position - org, hitpoint->position - org) > radius*radius)
			{
				return false;
			}
			if (dot(hitpoint->position - org, hitpoint->position - org) <  radius_inner * radius_inner)
			{
				return false;
			}


			if (!BlackHole)
			{
				double alp = 1.0;
				int xPos, yPos;
				if (hitpoint->material.disk_map && ((Circle*)this)->material()->texture)
				{
					DiscMapping discMap(0.0, radius, ((Circle*)this)->material()->texture->W(), ((Circle*)this)->material()->texture->H());
					discMap.MapUV(radius, (hitpoint->position - org).length(), hitpoint->u, hitpoint->v, xPos, yPos);

					if (xPos < 0 || yPos < 0)
					{
						return false;
					}
					Rgb t = ((Circle*)this)->material()->texture->cell(yPos, xPos);

					hitpoint->material.color.x = (double)t.r / 255.0;
					hitpoint->material.color.y = (double)t.g / 255.0;
					hitpoint->material.color.z = (double)t.b / 255.0;
					hitpoint->material.color = hitpoint->material.emission*hitpoint->material.color;

					alp = (double)t.alp / 255.0;
				}
				if (alp == 0.0)
				{
					return false;
				}
				if (alp < 0.99)
				{
					if (alp <= rnd->next01()) return false;
				}
			}

			//ç~íÖâ~î’
			if (BlackHole && this->blackhole_disk >= 0)
			{
				if (hitpoint->material.disk_map)
				{
					int xPos, yPos;
					if (discTexture[blackhole_disk] && discTexture[blackhole_disk]->GetImage())
					{
						double y[3];
						Vector3d p = BlackHole->Coordinate_transformation(hitpoint->position);

						Cartesian c(p.x, p.y, p.z);
						//Spherical s = c.ToSpherical();
						Spherical s = c.ToBoyerLindquist(BlackHole->a);
						

						//éñè€ÇÃínïΩñ Ç…ìûíB
						if (s.r < BlackHole->Rhor)
						{
							return false;
						}

						if (s.r > radius)
						{
							return false;
						}

						double alp = 1.0;
						if (s.r >= BlackHole->Rmstable)
						{
							//discMap->Map(s.r > BlackHole->Rdisk ? BlackHole->Rdisk:s.r, s.th, s.ph, xPos, yPos);
							discMap[blackhole_disk]->MapUV(radius, s.r, hitpoint->u, hitpoint->v, xPos, yPos);

							if (xPos < 0 || yPos < 0)
							{
								return false;
							}
							Rgb t = discTexture[blackhole_disk]->cell(yPos, xPos);

							hitpoint->material.color.x = (double)t.r / 255.0;
							hitpoint->material.color.y = (double)t.g / 255.0;
							hitpoint->material.color.z = (double)t.b / 255.0;
							//hitpoint->material.emission = ((Circle*)this)->material()->ibl_texture_coef*hitpoint->material.color;
							hitpoint->material.emission = hitpoint->material.emission*hitpoint->material.color;

							//alp = hitpoint->material.color.length()/3.0;
							alp = (double)t.alp / 255.0;
							alp = std::min( alp, hitpoint->material.color.length()/3.0);
							//if ( s.r > BlackHole->Rdisk*0.7 )
							//{
							//	alp *= exp(-pow(s.r - BlackHole->Rdisk*0.7, 1.2));
							//}
							//if (alp < 0.001) alp = 0.0;
						}
						else
						{
							return false;
						}


						if (alp == 0.0)
						{
							return false;
						}
						if (alp < 0.99)
						{
							if (alp <= rnd->next01()) return false;
						}
					}
				}
			}
			return true;
		}

		void MatrixTransformation(std::vector<Matrix4x4>& matrix)
		{
			base->MatrixTransformation(matrix);
			for (int j = 0; j < matrix.size(); j++)
			{
				org = matrix[j] * org;
				p[0] = matrix[j] * p[0];
				p[1] = matrix[j] * p[1];

				Matrix4x4 g = matrix[j].vector_coordinate_transformation();
				u_axis = g*u_axis;
				v_axis = g*v_axis;
				normal = g*normal;
			}
		}

	};
};

#endif
