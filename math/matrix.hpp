#ifndef MATH_MATRIX_HDR
#define MATH_MATRIX_HDR

#include <cmath>
#include <cassert>
#include <cstddef>
#include <algorithm>

#include "./math_defs.hpp"
#include "./tuple.hpp"

namespace newtonics {
	// forward-declarations
	namespace lib_math {
		template<typename type> struct t_point;
		template<typename type> struct t_vector;

		template<typename type> struct t_matrix {
		typedef t_point<type> m_point_type;
		typedef t_vector<type> m_vector_type;
		public:
			t_matrix<type>(const type* values) { set_values(values); }
			t_matrix<type>(const t_matrix<type>& m) { *this = m; }
			t_matrix<type>(
				const m_vector_type& x_vec = m_vector_type::x_axis_vector(),
				const m_vector_type& y_vec = m_vector_type::y_axis_vector(),
				const m_vector_type& z_vec = m_vector_type::z_axis_vector(),
				const m_vector_type& t_vec = m_vector_type::w_axis_vector()
			) {
				set_x_vector(x_vec);
				set_y_vector(y_vec);
				set_z_vector(z_vec);
				set_t_vector(t_vec);
			}

			t_matrix<type>& operator = (const t_matrix<type>& m) {
				set_values(m.get_values()); return *this;
			}


			bool operator == (const t_matrix<type>& m) const {
				unsigned int r = 0;

				for (unsigned int n = 0; n < MATH_MATRIX_SIZE; n += 4) {
					r += fp_eq(m_values[n + 0], m[n + 0], t_tuple<type>::eps_scalar());
					r += fp_eq(m_values[n + 1], m[n + 1], t_tuple<type>::eps_scalar());
					r += fp_eq(m_values[n + 2], m[n + 2], t_tuple<type>::eps_scalar());
					r += fp_eq(m_values[n + 3], m[n + 3], t_tuple<type>::eps_scalar());
				}

				return (r == MATH_MATRIX_SIZE);
			}

			bool operator != (const t_matrix<type>& m) const {
				return (!((*this) == m));
			}


			// matrix * point = point; 3x3 rotation + 1x3 translation
			// homogeneous w-coordinate of a point is 1 --> T included
			m_point_type operator * (const m_point_type& p) const {
				m_point_type tp;
				tp.x() = (m_values[0] * p.x()) + (m_values[4] * p.y()) + (m_values[ 8] * p.z()) + (m_values[12] * p.w());
				tp.y() = (m_values[1] * p.x()) + (m_values[5] * p.y()) + (m_values[ 9] * p.z()) + (m_values[13] * p.w());
				tp.z() = (m_values[2] * p.x()) + (m_values[6] * p.y()) + (m_values[10] * p.z()) + (m_values[14] * p.w());
				tp.w() = p.w();
				return tp;
			}

			// matrix * vector; interpret vector as point on unit-sphere
			// homogeneous w-coordinate of a vector is 0 --> T excluded
			//
			// only valid if matrix is orthonormal and vector has unit-length
			// in that case INVERSE(TRANSPOSE(R33)) == R33 and the vector can
			// be rotated directly (without requiring re-normalization due to
			// being translated as a point would be)
			//
			// aka "transform_vector"
			m_vector_type operator * (const m_vector_type& v) const {
				m_vector_type tv;
				tv.x() = (m_values[0] * v.x()) + (m_values[4] * v.y()) + (m_values[ 8] * v.z()); // same as n.inner_product(row[0])
				tv.y() = (m_values[1] * v.x()) + (m_values[5] * v.y()) + (m_values[ 9] * v.z()); // same as n.inner_product(row[1])
				tv.z() = (m_values[2] * v.x()) + (m_values[6] * v.y()) + (m_values[10] * v.z()); // same as n.inner_product(row[2])
				tv.w() = v.w();
				return tv;
			}



			// matrix * matrix = matrix
			t_matrix<type> operator * (const t_matrix<type>& m) const {
				t_matrix<type> r;

				// rows 0-2 with column 0 of m ( 0- 3)
				r[ 0] = (m_values[0] * m[ 0]) + (m_values[4] * m[ 1]) + (m_values[ 8] * m[ 2]) + (m_values[12] * m[ 3]);
				r[ 1] = (m_values[1] * m[ 0]) + (m_values[5] * m[ 1]) + (m_values[ 9] * m[ 2]) + (m_values[13] * m[ 3]);
				r[ 2] = (m_values[2] * m[ 0]) + (m_values[6] * m[ 1]) + (m_values[10] * m[ 2]) + (m_values[14] * m[ 3]);

				// rows 0-2 with column 1 of m ( 4- 7)
				r[ 4] = (m_values[0] * m[ 4]) + (m_values[4] * m[ 5]) + (m_values[ 8] * m[ 6]) + (m_values[12] * m[ 7]);
				r[ 5] = (m_values[1] * m[ 4]) + (m_values[5] * m[ 5]) + (m_values[ 9] * m[ 6]) + (m_values[13] * m[ 7]);
				r[ 6] = (m_values[2] * m[ 4]) + (m_values[6] * m[ 5]) + (m_values[10] * m[ 6]) + (m_values[14] * m[ 7]);

				// rows 0-2 with column 2 of m ( 8-11)
				r[ 8] = (m_values[0] * m[ 8]) + (m_values[4] * m[ 9]) + (m_values[ 8] * m[10]) + (m_values[12] * m[11]);
				r[ 9] = (m_values[1] * m[ 8]) + (m_values[5] * m[ 9]) + (m_values[ 9] * m[10]) + (m_values[13] * m[11]);
				r[10] = (m_values[2] * m[ 8]) + (m_values[6] * m[ 9]) + (m_values[10] * m[10]) + (m_values[14] * m[11]);

				// rows 0-2 with column 3 of m (12-15)
				r[12] = (m_values[0] * m[12]) + (m_values[4] * m[13]) + (m_values[ 8] * m[14]) + (m_values[12] * m[15]);
				r[13] = (m_values[1] * m[12]) + (m_values[5] * m[13]) + (m_values[ 9] * m[14]) + (m_values[13] * m[15]);
				r[14] = (m_values[2] * m[12]) + (m_values[6] * m[13]) + (m_values[10] * m[14]) + (m_values[14] * m[15]);

				return r;
			}

			t_matrix<type> operator * (const type values[MATH_MATRIX_SIZE]) const {
				return ((*this) * t_matrix<type>(&values[0]));
			}


			// matrix + matrix = matrix
			t_matrix<type> operator + (const t_matrix<type>& m) const {
				t_matrix<type> m_r;

				for (unsigned int n = 0; n < MATH_MATRIX_SIZE; n += 4) {
					m_r[n + 0] = m_values[n + 0] + m[n + 0];
					m_r[n + 1] = m_values[n + 1] + m[n + 1];
					m_r[n + 2] = m_values[n + 2] + m[n + 2];
					m_r[n + 3] = m_values[n + 3] + m[n + 3];
				}

				return m_r;
			}
			// matrix - matrix = matrix
			t_matrix<type> operator - (const t_matrix<type>& m) const {
				return ((*this) + (m * type(-1)));
			}

			/*
			t_matrix<type>  operator +  (const t_matrix<type>& m) const { t_matrix<type> mr(false); mr.add_values(m_values); mr.add_values(m.get_values()); return    mr; }
			t_matrix<type>  operator *  (const t_matrix<type>& m) const { t_matrix<type> mr(false); mr.add_values(m_values); mr.mul_values(m.get_values()); return    mr; }
			t_matrix<type>& operator += (const t_matrix<type>& m) const {                                                       add_values(m.get_values()); return *this; }
			t_matrix<type>& operator *= (const t_matrix<type>& m) const {                                                       mul_values(m.get_values()); return *this; }
			*/

			// matrix {*,+} matrix = matrix
			// slower than direct in-place updates but avoids duplication
			t_matrix<type>& operator *= (const t_matrix<type>& m) { return ((*this) = (*this) * m); }
			t_matrix<type>& operator += (const t_matrix<type>& m) { return ((*this) = (*this) + m); }


			// matrix * scalar = matrix
			t_matrix<type> operator * (const type s) const {
				t_matrix<type> m_r;

				for (unsigned int n = 0; n < MATH_MATRIX_SIZE; n += 4) {
					m_r[n + 0] = m_values[n + 0] * s;
					m_r[n + 1] = m_values[n + 1] * s;
					m_r[n + 2] = m_values[n + 2] * s;
					m_r[n + 3] = m_values[n + 3] * s;
				}

				return m_r;
			}

			t_matrix<type>& operator *= (const type s) { return ((*this) = (*this) * s); }


			type  operator [] (unsigned int idx) const { return m_values[idx]; }
			type& operator [] (unsigned int idx)       { return m_values[idx]; }

			// <idx> should be one of {0,4,8,12}
			const type* get_raw_vector(unsigned int idx) const { return &m_values[idx]; }
				  type* get_raw_vector(unsigned int idx)       { return &m_values[idx]; }
			const type* get_values() const { return &m_values[0]; }
			      type* get_values()       { return &m_values[0]; }


			// column accessors (NOTE: the t-vector is really a point!)
			m_vector_type get_x_vector() const { m_vector_type v; v.x() = m_values[ 0]; v.y() = m_values[ 1]; v.z() = m_values[ 2]; v.w() = m_values[ 3]; return v; }
			m_vector_type get_y_vector() const { m_vector_type v; v.x() = m_values[ 4]; v.y() = m_values[ 5]; v.z() = m_values[ 6]; v.w() = m_values[ 7]; return v; }
			m_vector_type get_z_vector() const { m_vector_type v; v.x() = m_values[ 8]; v.y() = m_values[ 9]; v.z() = m_values[10]; v.w() = m_values[11]; return v; }
			m_vector_type get_t_vector() const { m_vector_type v; v.x() = m_values[12]; v.y() = m_values[13]; v.z() = m_values[14]; v.w() = m_values[15]; return v; }

			// random accessors (<idx> should be one of {0,1,2,3})
			m_vector_type get_row_vector(unsigned int idx) const { return (m_vector_type(m_values[idx], m_values[idx + 4], m_values[idx + 8], m_values[idx + 12])); }
			m_vector_type get_col_vector(unsigned int idx) const { return (m_vector_type(m_values[(idx * 4) + 0], m_values[(idx * 4) + 1], m_values[(idx * 4) + 2], m_values[(idx * 4) + 3])); }


			t_matrix<type> translate(const m_vector_type& v) const {
				t_matrix<type> r = *this;
				r[12] += ((v.x() * r[0]) + (v.y() * r[4]) + (v.z() * r[ 8])); // same as tv.inner_product(rows[0])
				r[13] += ((v.x() * r[1]) + (v.y() * r[5]) + (v.z() * r[ 9])); // same as tv.inner_product(rows[1])
				r[14] += ((v.x() * r[2]) + (v.y() * r[6]) + (v.z() * r[10])); // same as tv.inner_product(rows[2])
				r[15] += ((v.x() * r[3]) + (v.y() * r[7]) + (v.z() * r[11])); // same as tv.inner_product(rows[3])
				return r;
			}
			t_matrix<type> scale(const m_vector_type& sv) const {
				t_matrix r = *this;

				r[ 0] *= sv.x(); r[ 4] *= sv.y();
				r[ 1] *= sv.x(); r[ 5] *= sv.y();
				r[ 2] *= sv.x(); r[ 6] *= sv.y();

				r[ 8] *= sv.z(); // r[12] *= sv.w();
				r[ 9] *= sv.z(); // r[13] *= sv.w();
				r[10] *= sv.z(); // r[14] *= sv.w();

				return r;
			}

			t_matrix<type> transpose_rotation() const {
				t_matrix<type> r = *this;
				std::swap(r[1], r[4]);
				std::swap(r[2], r[8]);
				std::swap(r[6], r[9]);
				return r;
			}
			t_matrix<type> transpose_translation() const {
				t_matrix<type> r = *this;
				std::swap(r[ 3], r[12]);
				std::swap(r[ 7], r[13]);
				std::swap(r[11], r[14]);
				return r;
			}
			t_matrix<type> transpose() const {
				t_matrix<type> r = *this;
				r.transpose_rotation_ref();
				r.transpose_translation_ref();
				return r;
			}

			// "affine" assumes matrix only does translation and rotation
			t_matrix<type> invert_affine() const {
				t_matrix<type> r = *this;

				// transpose the rotation
				r.transpose_rotation_ref();

				// get the inverse translation
				const m_vector_type t_pre = -r.get_t_vector();
				const m_vector_type t_inv = r * t_pre;

				// do the positional inversion
				r.set_t_vector(t_inv);
				return r;
			}
			// generalized inverse for non-orthonormal 4x4 matrices
			// A^-1 = (1 / det(A)) (C^T)_{ij} = (1 / det(A)) C_{ji}
			// where C is the matrix of cofactors
			//
			t_matrix<type> invert_general(const type eps = t_tuple<type>::eps_scalar()) const;

			t_matrix& translate_ref(const m_vector_type& tv) { return ((*this) = translate(tv)); }
			t_matrix& scale_ref(const m_vector_type& sv) { return ((*this) = scale(sv)); }

			t_matrix<type>& transpose_ref() { return ((*this) = transpose()); }
			t_matrix<type>& transpose_rotation_ref() { return ((*this) = transpose_rotation()); }
			t_matrix<type>& transpose_translation_ref() { return ((*this) = transpose_translation()); }
			t_matrix<type>& invert_affine_ref() { return ((*this) = invert_affine()); }



			t_matrix<type>& rotate_x_ref(type angle) { return ((*this) = rotate_x(angle)); }
			t_matrix<type>& rotate_y_ref(type angle) { return ((*this) = rotate_y(angle)); }
			t_matrix<type>& rotate_z_ref(type angle) { return ((*this) = rotate_z(angle)); }

			// these perform rotations around the (fixed aka external aka extrinsic) global coordinate axes
			// for any external rotation, multiply by the matrix R(C)*R(B)*R(A) if rotation order is A,B,C
			t_matrix<type>& rotate_xyz_ext_ref(const m_vector_type& angles) { return ((*this) = rotate_xyz_ext(angles)); }
			t_matrix<type>& rotate_yxz_ext_ref(const m_vector_type& angles) { return ((*this) = rotate_yxz_ext(angles)); }
			t_matrix<type>& rotate_zxy_ext_ref(const m_vector_type& angles) { return ((*this) = rotate_zxy_ext(angles)); }
			t_matrix<type>& rotate_zyx_ext_ref(const m_vector_type& angles) { return ((*this) = rotate_zyx_ext(angles)); }
			// these perform rotations around the (moving aka internal aka intrinsic) local coordinate axes
			// for any internal rotation, multiply by the matrix R(A)*R(B)*R(C) if rotation order is A,B,C
			// (or equivalently by the external matrix corresponding to the rotation order C,B,A)
			t_matrix<type>& rotate_xyz_int_ref(const m_vector_type& angles) { return ((*this) = rotate_xyz_int(angles)); }
			t_matrix<type>& rotate_yxz_int_ref(const m_vector_type& angles) { return ((*this) = rotate_yxz_int(angles)); }
			t_matrix<type>& rotate_zxy_int_ref(const m_vector_type& angles) { return ((*this) = rotate_zxy_int(angles)); }
			t_matrix<type>& rotate_zyx_int_ref(const m_vector_type& angles) { return ((*this) = rotate_zyx_int(angles)); }


			// rotate in X,Y,Z order [Rext = R(Z)*R(Y)*R(X), Rint = R(X)*R(Y)*R(Z)]
			t_matrix<type> rotate_xyz_ext(const m_vector_type angles) const {
				const t_matrix<type> m  = *this;
				const t_matrix<type> rx = std::move(compose_x_rotation(angles[P_AXIS_IDX]));
				const t_matrix<type> ry = std::move(compose_y_rotation(angles[Y_AXIS_IDX]));
				const t_matrix<type> rz = std::move(compose_z_rotation(angles[R_AXIS_IDX]));
				return (m * (rz * ry * rx));
			}
			t_matrix<type> rotate_xyz_int(const m_vector_type angles) const {
				return (rotate_zyx_ext(angles));
			}


			// rotate in Y,X,Z order [Rext = R(Z)*R(X)*R(Y), Rint = R(Y)*R(X)*R(Z)]
			t_matrix<type> rotate_yxz_ext(const m_vector_type angles) const {
				const t_matrix<type> m  = *this;
				const t_matrix<type> rx = std::move(compose_x_rotation(angles[P_AXIS_IDX]));
				const t_matrix<type> ry = std::move(compose_y_rotation(angles[Y_AXIS_IDX]));
				const t_matrix<type> rz = std::move(compose_z_rotation(angles[R_AXIS_IDX]));
				return (m * (rz * rx * ry));
			}
			t_matrix<type> rotate_yxz_int(const m_vector_type angles) const {
				return (rotate_zxy_ext(angles));
			}


			// rotate in Z,X,Y order [Rext = R(Y)*R(X)*R(Z), Rint = R(Z)*R(X)*R(Y)]
			t_matrix<type> rotate_zxy_ext(const m_vector_type angles) const {
				const t_matrix<type> m  = *this;
				const t_matrix<type> rx = std::move(compose_x_rotation(angles[P_AXIS_IDX]));
				const t_matrix<type> ry = std::move(compose_y_rotation(angles[Y_AXIS_IDX]));
				const t_matrix<type> rz = std::move(compose_z_rotation(angles[R_AXIS_IDX]));
				return (m * (ry * rx * rz));
			}
			t_matrix<type> rotate_zxy_int(const m_vector_type angles) const {
				return (rotate_yxz_ext(angles));
			}


			// rotate in Z,Y,X order [Rext = R(X)*R(Y)*R(Z), Rint = R(Z)*R(Y)*R(X)]
			t_matrix<type> rotate_zyx_ext(const m_vector_type angles) const {
				const t_matrix<type> m  = *this;
				const t_matrix<type> rx = std::move(compose_x_rotation(angles[P_AXIS_IDX]));
				const t_matrix<type> ry = std::move(compose_y_rotation(angles[Y_AXIS_IDX]));
				const t_matrix<type> rz = std::move(compose_z_rotation(angles[R_AXIS_IDX]));
				return (m * (rx * ry * rz));
			}
			t_matrix<type> rotate_zyx_int(const m_vector_type angles) const {
				return (rotate_xyz_ext(angles));
			}


			t_matrix<type> rotate_x(type angle) const { return ((*this) * compose_x_rotation(angle)); }
			t_matrix<type> rotate_y(type angle) const { return ((*this) * compose_y_rotation(angle)); }
			t_matrix<type> rotate_z(type angle) const { return ((*this) * compose_z_rotation(angle)); }


			// rotate in local or global YZ-plane; phi or alpha radians
			// constructs but does *not* apply Rx, so it can be chained
			static t_matrix<type> compose_x_rotation(type angle) {
				angle = clamp_angle_rad(angle);

				// Rx = ([X = [1, 0, 0], Y = [0, ca, -sa], Z = [0, sa, ca]])
				const type ca = std::cos(angle);
				const type sa = std::sin(angle);

				t_matrix<type> rm;
				rm[ 5] = +ca; rm[ 9] = +sa;
				rm[ 6] = -sa; rm[10] = +ca;
				return rm;
			}

			// rotate in local or global XZ-plane; theta or beta radians
			// constructs but does *not* apply Ry, so it can be chained
			static t_matrix<type> compose_y_rotation(type angle) {
				angle = clamp_angle_rad(angle);

				// Ry = ([X = [ca, 0, sa], Y = [0, 1, 0], Z = [-sa, 0, ca]])
				const type ca = std::cos(angle);
				const type sa = std::sin(angle);

				t_matrix<type> rm;
				rm[ 0] = +ca; rm[ 8] = -sa;
				rm[ 2] = +sa; rm[10] = +ca;
				return rm;
			}

			// rotate in local or global XY-plane; psi or gamma radians
			// constructs but does *not* apply Rz, so it can be chained
			static t_matrix<type> compose_z_rotation(type angle) {
				angle = clamp_angle_rad(angle);

				// Rz = ([X = [ca, -sa, 0], Y = [sa, ca, 0], Z = [0, 0, 1]])
				const type ca = std::cos(angle);
				const type sa = std::sin(angle);

				t_matrix<type> rm;
				rm[0] = +ca; rm[4] = +sa;
				rm[1] = -sa; rm[5] = +ca;
				return rm;
			}



			// rotation about any arbitrary axis, ie. in the plane crossing
			// the origin whose normal is given by <axis.x, axis.y, axis.z>
			// any such rotation can be decomposed into rotations about the
			// three component axes (by aligning the arbitrary axis with one
			// of the coordinate axes [here z], then rotating, then undoing
			// step 1)
			//
			// NOTE: rot_axis equal to Y is not supported, use rotate_y_ref
			t_matrix<type>& rotate_axis_ref(const m_vector_type& rot_axis, type angle) {
				return ((*this) = rotate_axis(rot_axis, angle));
			}

			t_matrix<type> rotate_axis(const m_vector_type& rot_axis, type angle) const {
				const m_vector_type vx = (rot_axis.outer_product(t_vector<type>::y_axis_vector())).normalize();
				const m_vector_type vy = (rot_axis.outer_product(                           vx  )).normalize();

				const t_matrix<type> m_fwd = std::move(t_matrix<type>(vx, vy, rot_axis));
				const t_matrix<type> m_inv = std::move(m_fwd.invert_affine());
				const t_matrix<type> m_rot = std::move(t_matrix<type>::compose_z_rotation(angle));

				return ((*this) * (m_fwd * m_rot * m_inv));
			}

			/*
			t_matrix<type>& rotate_axis_ref_tmp(const m_vector_type& rot_axis, type angle) {
				const type ca = std::cos(angle);
				const type sa = std::sin(angle);

				const m_vector_type x_axis = get_x_vector();
				const m_vector_type y_axis = get_y_vector();
				const m_vector_type z_axis = get_z_vector();

				// project rotation axis onto each principal axis
				const m_vector_type fwd_axis_x = rot_axis * x_axis.inner_product(rot_axis);
				const m_vector_type fwd_axis_y = rot_axis * y_axis.inner_product(rot_axis);
				const m_vector_type fwd_axis_z = rot_axis * z_axis.inner_product(rot_axis);

				// NOTE: rgt and up define the rotational plane
				const m_vector_type rgt_axis_x = x_axis - fwd_axis_x;
				const m_vector_type rgt_axis_y = y_axis - fwd_axis_y;
				const m_vector_type rgt_axis_z = z_axis - fwd_axis_z;

				// does not preserve orthonormality for non-principal rotation axes
				const m_vector_type up_axis_x = (rot_axis.outer_product(rgt_axis_x)).normalize();
				const m_vector_type up_axis_y = (rot_axis.outer_product(rgt_axis_y)).normalize();
				const m_vector_type up_axis_z = (rot_axis.outer_product(rgt_axis_z)).normalize();

				set_x_vector((fwd_axis_x + (rgt_axis_x * ca + up_axis_x * sa)));
				set_y_vector((fwd_axis_y + (rgt_axis_y * ca + up_axis_y * sa)));
				set_z_vector((fwd_axis_z + (rgt_axis_z * ca + up_axis_z * sa)));

				return *this;
			}
			*/



			// note: v.w() is ignored
			void set_x_vector(const m_vector_type& v) { m_values[ 0] = v.x(); m_values[ 1] = v.y(); m_values[ 2] = v.z(); m_values[ 3] = type(0); }
			void set_y_vector(const m_vector_type& v) { m_values[ 4] = v.x(); m_values[ 5] = v.y(); m_values[ 6] = v.z(); m_values[ 7] = type(0); }
			void set_z_vector(const m_vector_type& v) { m_values[ 8] = v.x(); m_values[ 9] = v.y(); m_values[10] = v.z(); m_values[11] = type(0); }
			void set_t_vector(const m_vector_type& v) { m_values[12] = v.x(); m_values[13] = v.y(); m_values[14] = v.z(); m_values[15] = type(1); }

			void set_row_vector(unsigned int idx, const m_vector_type& v) { m_values[(idx    ) + 0] = v.x(); m_values[(idx    ) + 4] = v.y(); m_values[(idx    ) + 8] = v.z(); m_values[(idx    ) + 12] = v.w(); }
			void set_col_vector(unsigned int idx, const m_vector_type& v) { m_values[(idx * 4) + 0] = v.x(); m_values[(idx * 4) + 1] = v.y(); m_values[(idx * 4) + 2] = v.z(); m_values[(idx * 4) +  3] = v.w(); }


			t_matrix<type>& set_xz_vectors_from_y_vector(const m_vector_type& vy) {
				set_y_vector(vy);

				// if y=<0,1,0> and z=<0,0,1>, sets x=<1,0,0> and z=<0,0,1>
				const m_vector_type& cvz = get_z_vector();
				const m_vector_type   vx = (vy.outer_product(cvz)).normalize();
				const m_vector_type   vz = (vx.outer_product( vy)).normalize();

				// avoid degeneracy if abs(inner_product(vy, cvz)) ~= 1
				assert(std::fabs(vy.inner_product(cvz)) < (type(1) - t_tuple<type>::eps_scalar()));

				set_x_vector(vx);
				set_z_vector(vz);
				return *this;
			}


			t_matrix<type>& add_values(const float values[MATH_MATRIX_SIZE]) {
				for (unsigned int i = 0; i < MATH_MATRIX_SIZE; i += 4) {
					m_values[i + 0] += values[i + 0];
					m_values[i + 1] += values[i + 1];
					m_values[i + 2] += values[i + 2];
					m_values[i + 3] += values[i + 3];
				}
				return *this;
			}

			t_matrix<type>& set_values(const float values[MATH_MATRIX_SIZE]) {
				for (unsigned int i = 0; i < MATH_MATRIX_SIZE; i += 4) {
					m_values[i + 0] = values[i + 0];
					m_values[i + 1] = values[i + 1];
					m_values[i + 2] = values[i + 2];
					m_values[i + 3] = values[i + 3];
				}
				return *this;
			}

			void set_null_values() {
				for (unsigned int n = 0; n < MATH_MATRIX_SIZE; n += 4) {
					m_values[n + 0] = type(0);
					m_values[n + 1] = type(0);
					m_values[n + 2] = type(0);
					m_values[n + 3] = type(0);
				}
			}
			void set_unit_values() {
				for (unsigned int n = 0; n < MATH_MATRIX_SIZE; n += 4) {
					m_values[n + 0] = type(n ==  0);
					m_values[n + 1] = type(n ==  4);
					m_values[n + 2] = type(n ==  8);
					m_values[n + 3] = type(n == 12);
				}
			}
			void print() const;

			unsigned int is_identity(const type eps = t_tuple<type>::eps_scalar()) const {
				const m_vector_type& xv = get_x_vector();
				const m_vector_type& yv = get_y_vector();
				const m_vector_type& zv = get_z_vector();
				const m_vector_type& tv = get_t_vector();

				unsigned int n = 0;

				n += ((!xv.equals(m_vector_type::x_axis_vector(), m_vector_type::x_axis_vector() * eps)) * (1 << 0));
				n += ((!yv.equals(m_vector_type::y_axis_vector(), m_vector_type::y_axis_vector() * eps)) * (1 << 1));
				n += ((!zv.equals(m_vector_type::z_axis_vector(), m_vector_type::z_axis_vector() * eps)) * (1 << 2));
				n += ((!tv.equals(m_vector_type::w_axis_vector(), m_vector_type::w_axis_vector() * eps)) * (1 << 3));

				return n;
			}
			unsigned int is_orthonormal(const type eps = t_tuple<type>::eps_scalar()) const {
				const m_vector_type& xv = get_x_vector();
				const m_vector_type& yv = get_y_vector();
				const m_vector_type& zv = get_z_vector();

				unsigned int n = 0;

				// test orthogonality
				n += ((std::fabs(xv.inner_product(yv)) >= eps) * (1 << 0));
				n += ((std::fabs(yv.inner_product(zv)) >= eps) * (1 << 1));
				n += ((std::fabs(xv.inner_product(zv)) >= eps) * (1 << 2));
				// test normality
				n += ((std::fabs(type(1) - xv.sq_magnit()) >= eps) * (1 << 3));
				n += ((std::fabs(type(1) - yv.sq_magnit()) >= eps) * (1 << 4));
				n += ((std::fabs(type(1) - zv.sq_magnit()) >= eps) * (1 << 5));

				return n;
			}


			m_vector_type get_angles_rh(const type eps = t_tuple<type>::eps_scalar()) const;
			m_vector_type get_angles_lh(const type eps = t_tuple<type>::eps_scalar()) const;

			//	matrix_outer(v, w) =
			//	  (  0 -v3  v2)   (w1)
			//	  ( v3   0 -v1) * (w2)
			//	  (-v2  v1   0)   (w3)
			//
			// NOTE: useful to generalize this to higher dimensions
			static m_vector_type outer_product(const m_vector_type& v, const m_vector_type& w) {
				t_matrix<type> m;
				m.set_x_vector(m_vector_type(type(0),   v.z(),  -v.y()));
				m.set_y_vector(m_vector_type( -v.z(), type(0),   v.x()));
				m.set_z_vector(m_vector_type(  v.y(),  -v.x(), type(0)));
				return (m * w);
			}

		private:
			// stores the elements in column-major order
			// m_values[0...3] represents the x-axis, etc
			type m_values[MATH_MATRIX_SIZE];
		};


		typedef t_matrix< float> t_matrix44f;
		typedef t_matrix<double> t_matrix44d;

		// aliases for less typing
		typedef t_matrix44f t_mat44f;
		typedef t_matrix44d t_mat44d;
	};
};

#endif

