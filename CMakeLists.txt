cmake_minimum_required(VERSION 3.1) #Проверка версии CMake.

set(PROJECT_NAME 000_Lebedev_AllLab) #Создаем обычную локальную переменную с именем проекта
set(CMAKE_CXX_STANDARD 11)
project(${PROJECT_NAME}) # Название проекта

set(PROBLEM_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/_bin CACHE STRING "Output directory")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG   ${PROBLEM_OUTPUT_DIRECTORY})
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE ${PROBLEM_OUTPUT_DIRECTORY})

set(GlobalSearch GlobalSearch) #Переменная с именем библиотеки
set(pugixml pugixml) #Переменная с именем библиотеки
set(src src)
option(USE_MPI "Set OFF to disable MPI" ON)

set(CMAKE_CONFIGURATION_TYPES "Debug;Release" CACHE STRING "Configs" FORCE) #Создается переменная CMAKE_CONFIGURATION_TYPES со значением "Debug;Release" и меткой/разделом Configs, переменная сохраниется в специальном файле доступном из всех файлов,
if(NOT CMAKE_BUILD_TYPE) # проверка имеющихся у компилятора вариантов сборки
  set(CMAKE_BUILD_TYPE Release) #Если нет то устанавливаем Release
endif()

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/pugixml)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/SRC/)

if(WIN32)
    option(gtest_force_shared_crt "" TRUE)
endif(WIN32)

#include( "${CMAKE_CURRENT_SOURCE_DIR}/cmake/Configure.cmake" )

############################### GTest ###############################
#add_subdirectory("${CMAKE_CURRENT_SOURCE_DIR}/3rdparty/gtest")
#include_directories("${CMAKE_SOURCE_DIR}/3rdparty/gtest/googletest/include")
#if( WIN32 )
 #   include_directories("${CMAKE_SOURCE_DIR}/3rdparty/gtest-mpi/win")
#else( WIN32 )
    #include_directories("${CMAKE_SOURCE_DIR}/3rdparty/gtest-mpi/linux")
#endif( WIN32 )

################################ MPI ################################

option(USE_MPI OFF)
if( USE_MPI )
    find_package( MPI )
    if( MPI_FOUND )
        include_directories( ${MPI_INCLUDE_PATH} )
    else( MPI_FOUND )
        set( USE_MPI OFF )
    endif( MPI_FOUND )
endif( USE_MPI )

add_subdirectory(pugixml)
add_subdirectory(src)
add_subdirectory(GlobalSearch)
add_subdirectory(BeautifulFunction)
add_subdirectory(Rastrigin2)
add_subdirectory(GKLS)
add_subdirectory(Grishagin)
add_subdirectory(StronginC3)


#target_link_libraries(src ${EXAMIN_MPFR_LIB} ${EXAMIN_MPIR_LIB} ${DISLIN_LIB})


#-----------------------------------------------------------------------------------
#file(GLOB hdrs "*.h*") #Добавляем в переменную hdrs все файлы с расширением .h
#file(GLOB srcs "*.cpp") #Добавляем в переменную srcs все файлы с расширением .cpp

#include_directories(lib1)
#add_executable(${PROJECT_NAME} ${srcs} ${hdrs}) #Создаем проект с именем из переменной ${PROJECT_NAME}, в него добавляются файлы из переменных ${srcs} и ${hdrs}


 #Собираем проект в поддиректории lib1
#target_link_libraries(${PROJECT_NAME} ${MyLib1}) #Подключаем к проекту с именем ${PROJECT_NAME} библиотеку с именем ${MyLib1}