#ifndef TESTLIB_H
#define TESTLIB_H

#include <iostream>
#include <list>
#include <string>
#include <map>
#include <vector>
#include <utility>



std::map<std::string,std::list<float>> zwshed_initial_c(const int dx, const int dy, const int dz, float* affs);

std::map<std::string,std::vector<double>> merge_with_stats(int dx,int dy, int dz, uint32_t * gt, float * rgn_graph,
                                        int rgn_graph_len, uint32_t * seg_in, uint32_t*counts, int counts_len, int thresh);

std::map<std::string,std::vector<double>> merge_no_stats(int dimX, int dimY, int dimZ, float * rgn_graph,
                                        int rgn_graph_len, uint32_t * seg_in, uint32_t*counts, int counts_len, int thresh);

std::map<std::string,std::list<float>> zwshed_initial_c_arb(const int dx, const int dy, const int dz, const uint32_t*node1,
                                               const uint32_t*node2, const float*edgeWeight, const int n_edge);

std::map<std::string,std::vector<double>> merge_with_stats_arb(int dx,int dy, int dz, uint32_t * gt, float * rgn_graph,
                                        int rgn_graph_len, uint32_t * seg_in, uint32_t*counts, int counts_len, int thresh);

std::map<std::string,std::vector<double>> merge_no_stats_arb(int dx,int dy, int dz, float * rgn_graph,
                                        int rgn_graph_len, uint32_t * seg_in, uint32_t*counts, int counts_len, int thresh);

#endif