#---------------------------------------------------------------------------
#
#  mimmo
#
#  Copyright (C) 2015-2017 OPTIMAD engineering Srl
#
#  -------------------------------------------------------------------------
#  License
#  This file is part of mimmo.
#
#  mimmo is free software: you can redistribute it and/or modify it
#  under the terms of the GNU Lesser General Public License v3 (LGPL)
#  as published by the Free Software Foundation.
#
#  mimmo is distributed in the hope that it will be useful, but WITHOUT
#  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
#  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
#  License for more details.
#
#  You should have received a copy of the GNU Lesser General Public License
#  along with mimmo. If not, see <http://www.gnu.org/licenses/>.
#
#---------------------------------------------------------------------------*/

#Specify the version being used as well as the language
cmake_minimum_required(VERSION 2.8)

# # Examples can be compiled only if all modules are enabled
# set(MODULE_LIST "${MIMMO_MODULE_LIST}")
# if (NOT ENABLE_MPI)
#     list(REMOVE_ITEM MODULE_LIST "${MPI_ONLY_MODULE_LIST}")
# endif ()

# areModulesEnabled("${MODULE_LIST}" ALL_MODULES_ENABLED)
# if (ALL_MODULES_ENABLED)
	option(BUILD_LAURATESTS "Create workspace for LAURA" ON)
# else ()
# 	unset(BUILD_LAURATESTS CACHE)
# 	return()
# endif ()

##NOTE###########
# Specify examples list according to module used. for instance, enable an example 
#using module utils verifying that a module utils is loaded 
################


# Add a target to generate the examples
foreach (MODULE_NAME IN LISTS MIMMO_MODULE_LIST)
	isModuleEnabled(${MODULE_NAME} MODULE_ENABLED)
	if (MODULE_ENABLED)
		addModuleIncludeDirectories(${MODULE_NAME})
	endif()
endforeach ()

if(BUILD_LAURATESTS)
	# List of examples
	
	file (GLOB LAURACODE *.cpp)
	
	set(LAURATEST_LIST "")
    list(APPEND LAURATEST_LIST "test_example")
    list(APPEND LAURATEST_LIST "testmain")
    list(APPEND LAURATEST_LIST "CalcoloOrto")
    list(APPEND LAURATEST_LIST "ExampleQuaternion")
    list(APPEND LAURATEST_LIST "ExamplePropagation")
    list(APPEND LAURATEST_LIST "ExampleDifferentCOR")
    list(APPEND LAURATEST_LIST "Naca0")


    
    foreach(LAURATEST_NAME IN LISTS LAURATEST_LIST)
        list(REMOVE_ITEM LAURACODE "${CMAKE_CURRENT_SOURCE_DIR}/${LAURATEST_NAME}.cpp")
    endforeach()
    
#     if (BUILD_SHARED_LIBS)
#         set(LAURA_LIBRARY_TYPE "SHARED")
#         set(LAURA_PIC_FLAG TRUE)
#     else()
         set(LAURA_LIBRARY_TYPE "STATIC")
#         set(LAURA_PIC_FLAG FALSE)
#     endif()
    set(LAURA_LIBRARY "mimmoEXTLAURA")
    add_library(${LAURA_LIBRARY} ${LAURA_LIBRARY_TYPE} ${LAURACODE})
    
#     set_target_properties(${LAURA_LIBRARY} PROPERTIES POSITION_INDEPENDENT_CODE ${PIC_FLAG})
#     set_target_properties(${LAURA_LIBRARY} PROPERTIES VERSION "${MIMMO_VERSION}"
#                                                  SOVERSION  "${MIMMO_MAJOR_VERSION}")
    #Rules to build the tests
	foreach(LAURATEST_NAME IN LISTS LAURATEST_LIST)
		set(LAURATEST_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/${LAURATEST_NAME}.cpp")
		list(APPEND LAURATESTS_SOURCES ${LAURATEST_SOURCES})

		add_executable(${LAURATEST_NAME} "${LAURATEST_SOURCES}")
		target_link_libraries(${LAURATEST_NAME} ${LAURA_LIBRARY})
		target_link_libraries(${LAURATEST_NAME} ${MIMMO_LIBRARY})
		target_link_libraries(${LAURATEST_NAME} ${MIMMO_EXTERNAL_LIBRARIES})
		install(TARGETS ${LAURATEST_NAME} DESTINATION doc/${CMAKE_PROJECT_NAME}-${MIMMO_VERSION}/laurawork)
	endforeach()

	add_custom_target(laura-tests DEPENDS ${LAURATEST_LIST})
	add_custom_target(clean-laura-tests COMMAND ${CMAKE_MAKE_PROGRAM} clean WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

	set (LAURATESTS_SOURCES ${LAURATESTS_SOURCES} PARENT_SCOPE)

    # Data files
    add_custom_command(
        TARGET "test_example" PRE_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/../geodata/planarquads_mimmo.vtu" "${CMAKE_CURRENT_BINARY_DIR}/geodata/planarquads_mimmo.vtu"
        COMMAND ${CMAKE_COMMAND} -E copy_if_different "${CMAKE_CURRENT_SOURCE_DIR}/../geodata/naca0012sq.vtu" "${CMAKE_CURRENT_BINARY_DIR}/geodata/naca0012sq.vtu"
    )
    

endif()




