#ifndef TESTLIB_H
#define TESTLIB_H

#include <iostream>
#include <list>
#include <string>
#include <map>
#include <vector>
#include <utility>

std::map<std::string,std::list<float>> calc_region_graph(int dx, int dy, int dz, int dcons, float* affs);

std::map<std::string,std::vector<double>> oneThresh_with_stats(int dx,int dy, int dz, int dcons, uint32_t * gt, float * affs, float * rgn_graph,
                                        int rgn_graph_len, uint32_t * seg_in, uint32_t*counts, int counts_len, int thresh,int eval);

std::map<std::string,std::vector<double>> oneThresh(int dimX, int dimY, int dimZ, int dcons, float* affs, float * rgn_graph,
                                        int rgn_graph_len, uint32_t * seg_in, uint32_t*counts, int counts_len, int thresh,int eval);

#endif