#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>

#include "boost/property_tree/ptree.hpp"
#include "boost/property_tree/json_parser.hpp"

#include "stdio.h"
#include "stdlib.h"

#include <map>
#include <unordered_map>
#include <boost/functional/hash.hpp>
#include <stack>
#include <algorithm>

#include "Point.h"
#include "json.hpp"

std::vector<std::string> spliter(const std::string& str){
    std::vector<std::string> splited;
    std::string element;
    for (auto x : str){
        if (x == ',')
        {
            splited.push_back(element);
            element = "";
        }
        else {
            element = element + x;
        }
    }
    splited.push_back(element);
    return splited;
}

void read(const char *file_in, std::unordered_map<std::string, std::vector<Point>> &v,
          std::unordered_map<std::string, std::string>& constr,
          std::unordered_map<std::string, unsigned int>& floors ) {
    std::string line = "";
    std::ifstream file(file_in, std::ifstream::in);
    if (!file) { std::cerr << "Input file not found. Check the relative file path\n"; }
    int count = 0;
    while (!file.eof()) {
        std::getline(file, line);
        if (count == 0){count++; continue;}
        std::istringstream iss(line);
        std::vector<std::string> splited_line = spliter(line);
        std::vector<Point> verts;
        for (unsigned int i = 16; i<splited_line.size(); i=i+2 ) {

            size_t j = 0;
            for ( ; j < splited_line[i].length(); j++ ){ if ( isdigit(splited_line[i][j]) ) break; }
            splited_line[i] = splited_line[i].substr(j, splited_line[i].length() - j );
            double id_x = atof(splited_line[i].c_str());

            j=0;
            for ( ; j < splited_line[i+1].length(); j++ ){ if ( isdigit(splited_line[i+1][j]) ) break; }
            splited_line[i+1] = splited_line[i+1].substr(j, splited_line[i+1].length() - j );
            double id_y = atof(splited_line[i+1].c_str());
            double z = stod(splited_line[13]);
            std::pair<double,double> p;
            p.first = id_x; p.second = id_y;
            double roofz = stod(splited_line[12]);


            verts.push_back(Point(id_x,id_y,z,roofz));
                //v.push_back(Point(std::stof(splited_line[1]), std::stof(splited_line[2]), std::stof(splited_line[3])));

        }
        if (line == ""){continue;}
        std::string id = splited_line[1];
        v[id] = verts;
        constr[id]=splited_line[8];
        floors[id]=floor(verts.back().z_r/3);
    }
}

std::vector<double> cornerpoints(std::vector<std::vector<double>> v, const std::string& minmax){
    if (minmax == "max"){
        float maxx = -2000.0;
        float maxy = -2000.0;
        float maxz = -2000.0;
        for (auto & i : v) {
            if (i[0] >= maxx) maxx = i[0];
            if (i[1] >= maxy) maxy = i[1];
            if (i[2] >= maxz) maxz = i[2];
        }
        return {maxx,maxy,maxz}; //+1 to go out of bounds
    }
    else if (minmax == "min") {
        float minx = 2000.0;
        float miny = 2000.0;
        float minz = 2000.0;
        for (auto & i : v) {
            if (i[0] <= minx) minx = i[0];
            if (i[1] <= miny) miny = i[1];
            if (i[2] <= minz) minz = i[2];
        }
        return {minx,miny,minz}; //-1 to go out of bounds
    }
}

void populate(std::unordered_map<std::string, std::unordered_map<unsigned int, std::vector<Point>>>& b,std::unordered_map<std::string, std::vector<Point>>& v){

    for (auto const& i : v){
        std::unordered_map<unsigned int, std::vector<Point>> second;
        for (unsigned int j = 0; j < i.second.size()+2; j++) { //faces
            std::vector<Point> vertices;
            if (j == 0) {
                for (unsigned int k = 0; k < i.second.size(); k++) {
                    vertices.push_back(i.second[k]);
                }
            }
            else if ( j == 1){
                    for (unsigned int k = 0; k < i.second.size(); k++) {

                        vertices.push_back(Point(i.second[k].x,i.second[k].y,i.second[k].z_r,i.second[k].z));
                    }
                }
            else{
                //if (j-2+1 == i.second.size()){break;}
                    //auto z = i.second[k].z + roof[std::make_pair(i.second[k].x,i.second[k].y)];
                    auto bottom1 = b[i.first].find(0)->second[j-2];
                vertices.push_back(bottom1);
                    auto bottom2 = b[i.first].find(0)->second[j-2+1];
                vertices.push_back(bottom2);
                    auto up1 = b[i.first].find(1)->second[j-2];
                vertices.push_back(up1);
                    auto up2 = b[i.first].find(1)->second[j-2+1];
                vertices.push_back(up2);
            }



            second[j] = vertices;
            b[i.first] = second;
        }
    }
}

int main(int argc, const char * argv[]) {
    const char *file_in = "/home/konstantinos/Desktop/TUDelft-Courses/Q3/GEO1004/hw3/data/tudcampus.csv";
    const char *file_out = "/home/konstantinos/Desktop/TUDelft-Courses/Q3/GEO1004/hw3/data/TU.json";
    std::unordered_map<std::string, std::vector<Point>> vertices;
    std::unordered_map<std::string, std::unordered_map<unsigned int, std::vector<Point>>> buildings;
    std::unordered_map<std::string, std::string> const_year;
    std::unordered_map<std::string, unsigned int> storeys;

    read(file_in,vertices,const_year, storeys);


    populate(buildings,vertices);
    auto a= storeys["0503100000020276"];
    int count=0;

  return 0;
}
