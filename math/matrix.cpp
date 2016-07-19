#include <cstdio>

#include "./matrix.hpp"
#include "./point.hpp"
#include "./vector.hpp"

namespace newtonics {
	namespace lib_math {
		#if 0
		template<> t_matrix44f t_matrix44f::invert_general(const float eps) const {
			t_matrix44f mat;
			t_matrix44f& cofac = mat;

			for (unsigned int i = 0; i < MATH_MATRIX_SIZE; i++) {
				for (unsigned int j = 0; j < MATH_MATRIX_SIZE; j++) {
					cofac.m_values[i][j] = calc_cofac(m_values, i, j);
				}
			}

			const float det =
				(m_values[0][0] * cofac.m_values[0][0]) +
				(m_values[0][1] * cofac.m_values[0][1]) +
				(m_values[0][2] * cofac.m_values[0][2]) +
				(m_values[0][3] * cofac.m_values[0][3]);

			if (det < eps) {
				// <this> is a singular matrix
				return mat;
			}

			mat *= (1.0f / det);
			mat.transpose_ref();
			return mat;
		}
		#endif


		// extracts the rotation angles from a matrix
		//
		// this assumes a right-handed coordinate system; if called
		// on a left-handed matrix the angles {a,b,c} returned here
		// will correspond to the effects of rotate_zyx(-a, -b, -c)
		// instead of rotate_xyz(a, b, c) in the RH case
		//
		// however, since
		//
		//   1)           rotate_zyx(-a, -b, -c)  == transpose(rotate_xyz(a, b, c))
		//   2) transpose(rotate_zyx(-a, -b, -c)) ==           rotate_xyz(a, b, c)
		//
		// and any internal rotation A(a)B(b)C(c) equals the external
		// rotation C(c)B(b)A(a) we can easily convert them back to a
		// left-handed form (left-handed coordinate system is used by
		// default)
		//
		// all angles are in radians and returned in PYR order
		//
		template<> t_vector4f t_matrix44f::get_angles_rh(const float eps) const {
			t_vector4f angles[2];

			float cos_yaw[2] = {0.0f, 0.0f};
			float ang_sum[2] = {0.0f, 0.0f};

			if (square(eps) > std::fabs(m_values[0 * 4 + 2] + 1.0f)) {
				// x.z == -1 (yaw=PI/2) means gimbal lock between X and Z
				angles[0][R_AXIS_IDX] = (       0.0f);
				angles[0][Y_AXIS_IDX] = (M_PI * 0.5f);
				angles[0][P_AXIS_IDX] = (angles[0][R_AXIS_IDX] + std::atan2(m_values[1 * 4 + 0], m_values[2 * 4 + 0]));

				return angles[0];
			}

			if (square(eps) > std::fabs(m_values[0 * 4 + 2] - 1.0f)) {
				// x.z == 1 (yaw=-PI/2) means gimbal lock between X and Z
				angles[0][R_AXIS_IDX] =  (       0.0f);
				angles[0][Y_AXIS_IDX] = -(M_PI * 0.5f);
				angles[0][P_AXIS_IDX] = (-angles[0][R_AXIS_IDX] + std::atan2(-m_values[1 * 4 + 0], -m_values[2 * 4 + 0]));

				return angles[0];
			}

			// yaw angles (theta)
			//
			//   angles[i][P] :=   psi := Pitch := X-angle
			//   angles[i][Y] := theta :=   Yaw := Y-angle
			//   angles[i][R] :=   phi :=  Roll := Z-angle
			//
			angles[0][Y_AXIS_IDX] = -std::asin(m_values[0 * 4 + 2]);
			angles[1][Y_AXIS_IDX] = (M_PI - angles[0][Y_AXIS_IDX]);

			// yaw cosines
			cos_yaw[0] = std::cos(angles[0][Y_AXIS_IDX]);
			cos_yaw[1] = std::cos(angles[1][Y_AXIS_IDX]);

			// psi angles (pitch)
			angles[0][P_AXIS_IDX] = std::atan2((m_values[1 * 4 + 2] / cos_yaw[0]), (m_values[2 * 4 + 2] / cos_yaw[0]));
			angles[1][P_AXIS_IDX] = std::atan2((m_values[1 * 4 + 2] / cos_yaw[1]), (m_values[2 * 4 + 2] / cos_yaw[1]));
			// phi angles (roll)
			angles[0][R_AXIS_IDX] = std::atan2((m_values[0 * 4 + 1] / cos_yaw[0]), (m_values[0 * 4 + 0] / cos_yaw[0]));
			angles[1][R_AXIS_IDX] = std::atan2((m_values[0 * 4 + 1] / cos_yaw[1]), (m_values[0 * 4 + 0] / cos_yaw[1]));

			ang_sum[0] = (m_vector_type::ones_vector()).inner_product(angles[0].abs()); // |p0|+|y0|+|r0|
			ang_sum[1] = (m_vector_type::ones_vector()).inner_product(angles[1].abs()); // |p1|+|y1|+|r1|

			// two solutions exist; choose the "shortest" rotation
			return angles[ang_sum[0] > ang_sum[1]];
		}

		template<> t_vector4f t_matrix44f::get_angles_lh(const float eps) const {
			t_matrix44f matrix = *this;
			t_vector4f angles;
			// we are given a left-handed matrix (either xyz_int or zyx_ext), but the above
			// function expects it to be right-handed which is solved simply by transposing
			// since (A(a)*B(b)*C(c))^T = (C(-c)^T)*(B(-b)^T)*(A(-a)^T)
			matrix = std::move(matrix.transpose());
			angles = std::move(matrix.get_angles_rh(eps));
			return angles;
		}


		template<> void t_matrix44f::print() const {
			printf("\n X      Y      Z      T\n");
			for (unsigned int n = 0; n < MATH_MATRIX_SIZE; n += 4) { printf(" %.3f ", m_values[n]); } printf("\n"); // 1st row
			for (unsigned int n = 1; n < MATH_MATRIX_SIZE; n += 4) { printf(" %.3f ", m_values[n]); } printf("\n"); // 2nd row
			for (unsigned int n = 2; n < MATH_MATRIX_SIZE; n += 4) { printf(" %.3f ", m_values[n]); } printf("\n"); // 3rd row
			for (unsigned int n = 3; n < MATH_MATRIX_SIZE; n += 4) { printf(" %.3f ", m_values[n]); } printf("\n"); // 4th row
			printf("\n");
		}
	};
};

