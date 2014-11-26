include /user/kmarcus2/Elemental/install/conf/elemvariables

kyle-coifman-elemental: kyle-coifman-elemental.cpp coifman.h libMPIdebug.h
	${CXX} ${ELEM_COMPILE_FLAGS} $< -o $@ ${ELEM_LINK_FLAGS} ${ELEM_LIBS}
