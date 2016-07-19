#include <cassert>

#include <GL/gl.h>
#include <GL/glu.h>

#include "opengl_wrapper.hpp"

t_gl_wrapper_state gl;

void t_gl_wrapper_state::InitClientArrays(const t_vertex* verts) {
	assert(      verts     != nullptr);
	assert(m_cur_vert_pntr != nullptr);

	// do not check if verts == m_cur_vert_pntr; is already true
	// on the first End() before any call to gl*Pointer has been
	// made
	m_prv_vert_pntr = m_cur_vert_pntr;
	m_cur_vert_pntr = verts;

	glVertexPointer  (  sizeof(verts->pxyzw) / sizeof(  float),   GL_FLOAT        , sizeof(t_vertex), verts->pxyzw);
	glNormalPointer  (/*sizeof(verts->nxyz ) / sizeof(  float),*/ GL_FLOAT        , sizeof(t_vertex), verts->nxyz );
	glTexCoordPointer(  sizeof(verts->stuv ) / sizeof(  float),   GL_FLOAT        , sizeof(t_vertex), verts->stuv );
	glColorPointer   (  sizeof(verts->rgba ) / sizeof(uint8_t),   GL_UNSIGNED_BYTE, sizeof(t_vertex), verts->rgba );
}

void t_gl_wrapper_state::EnableClientArrays() const {
	glEnableClientState(GL_VERTEX_ARRAY);
	glEnableClientState(GL_NORMAL_ARRAY);
	glEnableClientState(GL_TEXTURE_COORD_ARRAY);
	glEnableClientState(GL_COLOR_ARRAY);
}

void t_gl_wrapper_state::DisableClientArrays() const {
	glDisableClientState(GL_VERTEX_ARRAY);
	glDisableClientState(GL_NORMAL_ARRAY);
	glDisableClientState(GL_TEXTURE_COORD_ARRAY);
	glDisableClientState(GL_COLOR_ARRAY);
}

void t_gl_wrapper_state::DrawElements(unsigned int mode, unsigned int count, const t_vertex* verts, const t_index* inds) {
	InitClientArrays(verts);
	glDrawElements(mode, count, GL_UNSIGNED_INT, inds);
	InitClientArrays(GetArrayPointer()); // reset
}

void t_gl_wrapper_state::End() {
	// too risky, just call ICA
	// assert(m_cur_vert_pntr == &m_verts[0]);

	InitClientArrays(&m_verts[0]);
	glDrawArrays(m_cur_prim_type, 0, m_cur_vert_indx);
	InitClientArrays(m_prv_vert_pntr);
}

