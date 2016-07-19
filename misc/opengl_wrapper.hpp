#ifndef NEWTONICS_OPENGL_WRAPPER_HDR
#define NEWTONICS_OPENGL_WRAPPER_HDR

#include <cstdint>
// #include <array>
#include <vector>

// wrapper state for immediate-mode functions; should
// never be used except to lift old GL code out of IM
// with minimal effort
//
// template<unsigned int va_size = 16384>
struct t_gl_wrapper_state {
public:
	struct t_vertex {
		float   pxyzw[4]; // position
		float   nxyz [3]; // normal
		float   stuv [4]; // texcoors
		uint8_t rgba [4]; // color
	};

	typedef unsigned int t_index;

public:
	t_gl_wrapper_state(unsigned int va_size = 16384) {
		m_cur_vert_indx = 0;
		m_cur_prim_type = 0;
		m_max_num_verts = va_size;

		// assume each vertex is referenced by three triangles
		m_verts.resize(m_max_num_verts + 1);
		m_inds.resize(m_max_num_verts * 3, 0);

		for (t_vertex& vert: m_verts) {
			vert.pxyzw[0] = 0.0f;
			vert.pxyzw[1] = 0.0f;
			vert.pxyzw[2] = 0.0f;
			vert.pxyzw[3] = 0.0f;

			vert.nxyz[0] = 0.0f;
			vert.nxyz[1] = 0.0f;
			vert.nxyz[2] = 0.0f;

			vert.stuv[0] = 0.0f;
			vert.stuv[1] = 0.0f;
			vert.stuv[2] = 0.0f;
			vert.stuv[3] = 0.0f;

			vert.rgba[0] = 0;
			vert.rgba[1] = 0;
			vert.rgba[2] = 0;
			vert.rgba[3] = 0;
		}

		// prv is set on first call to ICA
		m_cur_vert_pntr = GetArrayPointer();
		m_prv_vert_pntr = nullptr;
	}

	// convenience wrapper; barely justified
	void DrawElements(unsigned int mode, unsigned int count, const t_vertex* verts, const t_index* inds);
	// only has to be called at startup or in End() with GetArrayPointer() as arg
	// (if indexed drawing is also performed, this is taken care of automatically)
	void InitClientArrays(const t_vertex* verts);

	void EnableClientArrays() const;
	void DisableClientArrays() const;

	void Begin(int prim_type) {
		m_cur_vert_indx = 0;
		m_cur_prim_type = prim_type;
	}

	void End();


	// TODO: MultiTexCoord*{f,d,i,s}
	void TexCoord1f(float s                           ) { TexCoord4f(s, 0.0f, 0.0f, 0.0f); }
	void TexCoord2f(float s, float t                  ) { TexCoord4f(s,    t, 0.0f, 0.0f); }
	void TexCoord3f(float s, float t, float u         ) { TexCoord4f(s,    t,    u, 0.0f); }
	void TexCoord4f(float s, float t, float u, float v) {
		t_vertex& vert = m_verts[m_max_num_verts];
		vert.stuv[0] = s;
		vert.stuv[1] = t;
		vert.stuv[2] = u;
		vert.stuv[3] = v;
	}

	// TODO: SecondaryColor*{f,d,[u]b,[u]s,[u]i}
	void Color1f(float r                           ) { Color4ub(r * 255.0f,       0.0f,       0.0f,     255.0f); }
	void Color2f(float r, float g                  ) { Color4ub(r * 255.0f, g * 255.0f,       0.0f,     255.0f); }
	void Color3f(float r, float g, float b         ) { Color4ub(r * 255.0f, g * 255.0f, b * 255.0f,     255.0f); }
	void Color4f(float r, float g, float b, float a) { Color4ub(r * 255.0f, g * 255.0f, b * 255.0f, a * 255.0f); }

	void Color1ub(uint8_t r                                 ) { Color4ub(r, 0, 0, 255); }
	void Color2ub(uint8_t r, uint8_t g                      ) { Color4ub(r, g, 0, 255); }
	void Color3ub(uint8_t r, uint8_t g, uint8_t b           ) { Color4ub(r, g, b, 255); }
	void Color4ub(uint8_t r, uint8_t g, uint8_t b, uint8_t a) {
		t_vertex& vert = m_verts[m_max_num_verts];
		vert.rgba[0] = r;
		vert.rgba[1] = g;
		vert.rgba[2] = b;
		vert.rgba[3] = a;
	}

	// TODO: Normal{b,d,i,s}
	void Normal3f(float x, float y, float z) {
		t_vertex& vert = m_verts[m_max_num_verts];
		vert.nxyz[0] = x;
		vert.nxyz[1] = y;
		vert.nxyz[2] = z;
	}

	// TODO: Vertex{b,d,i,s}
	void Vertex2f(float x, float y                  ) { Vertex4f(x, y, 0.0f, 1.0f); }
	void Vertex3f(float x, float y, float z         ) { Vertex4f(x, y,    z, 1.0f); }
	void Vertex4f(float x, float y, float z, float w) {
		t_vertex& vert = m_verts[m_max_num_verts];
		vert.pxyzw[0] = x;
		vert.pxyzw[1] = y;
		vert.pxyzw[2] = z;
		vert.pxyzw[3] = w;

		// copy attributes, advance index but never overflow ('&' needs a POT-1)
		m_verts[(m_cur_vert_indx++) % m_max_num_verts] = m_verts[m_max_num_verts];
	}

	const t_vertex* GetArrayPointer() const { return &m_verts[0]; }
	      t_vertex* GetArrayPointer()       { return &m_verts[0]; }

	const t_index* GetIndexPointer() const { return &m_inds[0]; }
	      t_index* GetIndexPointer()       { return &m_inds[0]; }

private:
	// last vertex holds current attributes
	// std::array<t_vertex, va_size + 1> m_verts;
	std::vector<t_vertex> m_verts;
	std::vector<t_index> m_inds;

	const t_vertex* m_cur_vert_pntr;
	const t_vertex* m_prv_vert_pntr;

	unsigned int m_max_num_verts;
	unsigned int m_cur_vert_indx;
	unsigned int m_cur_prim_type;
};

extern t_gl_wrapper_state gl;

#ifndef GL_VERSION
#error "wrapper must be included after <GL/gl.h>"
#endif

#undef glBegin
#undef glEnd
#undef glVertex2f
#undef glVertex3f
#undef glVertex4f
#undef glNormal3f
#undef glTexCoord1f
#undef glTexCoord2f
#undef glTexCoord3f
#undef glTexCoord4f
#undef glColor1f
#undef glColor2f
#undef glColor3f
#undef glColor4f
#undef glColor1ub
#undef glColor2ub
#undef glColor3ub
#undef glColor4ub

// re-route all IM calls; also causes nice compilation
// error if gl.h is ever included after us by accident
#define glBegin gl.Begin
#define glEnd gl.End
#define glVertex2f gl.Vertex2f
#define glVertex3f gl.Vertex3f
#define glVertex4f gl.Vertex4f
#define glNormal3f gl.Normal3f
#define glTexCoord1f gl.TexCoord1f
#define glTexCoord2f gl.TexCoord2f
#define glTexCoord3f gl.TexCoord3f
#define glTexCoord4f gl.TexCoord4f
#define glColor1f gl.Color1f
#define glColor2f gl.Color2f
#define glColor3f gl.Color3f
#define glColor4f gl.Color4f
#define glColor1ub gl.Color1ub
#define glColor2ub gl.Color2ub
#define glColor3ub gl.Color3ub
#define glColor4ub gl.Color4ub

#endif

