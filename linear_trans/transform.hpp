#ifndef KMUCS_GRAPHICS_TRANSFORM_HPP
#define KMUCS_GRAPHICS_TRANSFORM_HPP

#include <cmath>
#include "vec.hpp"
#include "mat.hpp"
#include "operator.hpp"

namespace kmuvcl
{
	namespace math
	{
#ifndef M_PI
		const float M_PI = 3.14159265358979323846f;
#endif

		template <typename T>
		mat<4, 4, T> translate(T dx, T dy, T dz)
		{
			mat<4, 4, T> translateMat;

			// TODO: Fill up this function properly 
			for (int i = 0; i < 4; ++i)
			{
				translateMat(i, i) = 1;
			}
			translateMat(0, 3) = dx;
			translateMat(1, 3) = dy;
			translateMat(2, 3) = dz;

			return translateMat;
		}

		template <typename T>
		mat<4, 4, T> rotate(T angle, T x, T y, T z)
		{
			mat<4, 4, T> rotateMat;

			// TODO: Fill up this function properly 
			T radian = (angle * M_PI) / 180;
			T cos_radian = cos(radian);
			T sin_radian = sin(radian);
			T scalar = sqrt(pow(x, 2) + pow(y, 2) + pow(z, 2));
			T ux = x / scalar;
			T uy = y / scalar;
			T uz = z / scalar;

			rotateMat(0, 0) = cos_radian + pow(ux, 2)*(1 - cos_radian);
			rotateMat(0, 1) = ux * uy*(1 - cos_radian) - uz * sin_radian;
			rotateMat(0, 2) = ux * uz*(1 - cos_radian) + uy * sin_radian;

			rotateMat(1, 0) = ux * uy*(1 - cos_radian) + uz * sin_radian;
			rotateMat(1, 1) = cos_radian + pow(uy, 2)*(1 - cos_radian);
			rotateMat(1, 2) = uy * uz*(1 - cos_radian) - ux * sin_radian;

			rotateMat(2, 0) = ux * uz*(1 - cos_radian) - uy * sin_radian;
			rotateMat(2, 1) = uy * uz*(1 - cos_radian) + ux * sin_radian;
			rotateMat(2, 2) = cos_radian + pow(uz, 2)*(1 - cos_radian);

			rotateMat(3, 3) = 1;
			return rotateMat;
		}

		template<typename T>
		mat<4, 4, T> scale(T sx, T sy, T sz)
		{
			mat<4, 4, T> scaleMat;

			// TODO: Fill up this function properly
			scaleMat(0, 0) = sx;
			scaleMat(1, 1) = sy;
			scaleMat(2, 2) = sz;
			scaleMat(3, 3) = 1;

			return scaleMat;
		}

		template<typename T>
		mat<4, 4, T> lookAt(T eyeX, T eyeY, T eyeZ, T centerX, T centerY, T centerZ, T upX, T upY, T upZ)
		{
			mat<4, 4, T> viewMat;

			// TODO: Fill up this function properly 
			
			T cam_pos_x = eyeX;
			T cam_pos_y = eyeY;
			T cam_pos_z = eyeZ;

			T scalar = sqrt(pow((eyeX - centerX), 2) + pow((eyeY - centerY), 2) + pow((eyeZ - centerZ), 2));
			T cam_z_axis_x = (eyeX - centerX) / scalar;
			T cam_z_axis_y = (eyeY - centerY) / scalar;
			T cam_z_axis_z = (eyeZ - centerZ) / scalar;

			T ux = upY * cam_z_axis_z - upZ * cam_z_axis_y;
			T uy = upZ * cam_z_axis_x - upX * cam_z_axis_z;
			T uz = upX * cam_z_axis_y - upY * cam_z_axis_x;

			scalar = sqrt(pow(ux, 2) + pow(uy, 2) + pow(uz, 2));
			T cam_x_axis_x = ux/scalar;
			T cam_x_axis_y = uy/scalar;
			T cam_x_axis_z = uz/scalar;

			ux = cam_z_axis_y * cam_x_axis_z - cam_z_axis_z * cam_x_axis_y;
			uy = cam_z_axis_z * cam_x_axis_x - cam_z_axis_x * cam_x_axis_z;
			uz = cam_z_axis_x * cam_x_axis_y - cam_z_axis_y * cam_x_axis_x;
			scalar = sqrt(pow(ux, 2) + pow(uy, 2) + pow(uz, 2));

			T cam_y_axis_x = ux / scalar;
			T cam_y_axis_y = uy / scalar;
			T cam_y_axis_z = uz / scalar;

			mat<4, 4, T> viewMatLeft;
			mat<4, 4, T> viewMatRight;
			
			viewMatLeft(0, 0) = cam_x_axis_x;
			viewMatLeft(0, 1) = cam_x_axis_y;
			viewMatLeft(0, 2) = cam_x_axis_z;
			viewMatLeft(1, 0) = cam_y_axis_x;
			viewMatLeft(1, 1) = cam_y_axis_y;
			viewMatLeft(1, 2) = cam_y_axis_z;
			viewMatLeft(2, 0) = cam_z_axis_x;
			viewMatLeft(2, 1) = cam_z_axis_y;
			viewMatLeft(2, 2) = cam_z_axis_z;
			viewMatLeft(3, 3) = 1;

			for (int i = 0; i < 4; ++i)
			{
				viewMatRight(i, i) = 1;
			}
			viewMatRight(0, 3) = -cam_pos_x;
			viewMatRight(1, 3) = -cam_pos_y;
			viewMatRight(2, 3) = -cam_pos_z;

			/*std::cout << "viewMatLeft : " << std::endl;
			std::cout << viewMatLeft << std::endl;
			std::cout << "viewMatRight : " << std::endl;
			std::cout << viewMatRight << std::endl;*/

			viewMat = viewMatLeft * viewMatRight;

			return viewMat;
		}

		template<typename T>
		mat<4, 4, T> ortho(T left, T right, T bottom, T top, T nearVal, T farVal)
		{
			mat<4, 4, T> orthoMat;

			// TODO: Fill up this function properly 
			orthoMat(0, 0) = 2 / (right - left);
			orthoMat(1, 1) = 2 / (top - bottom);
			orthoMat(2, 2) = -2 / (farVal - nearVal);
			orthoMat(3, 3) = 1;
			orthoMat(0, 3) = -(right + left) / (right - left);
			orthoMat(1, 3) = -(top + bottom) / (top - bottom);
			orthoMat(2, 3) = -(farVal + nearVal) / (farVal - nearVal);

			return orthoMat;
		}

		template<typename T>
		mat<4, 4, T> frustum(T left, T right, T bottom, T top, T nearVal, T farVal)
		{
			mat<4, 4, T> frustumMat;

			// TODO: Fill up this function properly 

			return frustumMat;
		}

		template<typename T>
		mat<4, 4, T> perspective(T fovy, T aspect, T zNear, T zFar)
		{
			T  right = 0;
			T  top = 0;

			// TODO: Fill up this function properly 

			return frustum(-right, right, -top, top, zNear, zFar);
		}
	}
}
#endif
