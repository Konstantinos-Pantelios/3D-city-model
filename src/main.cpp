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
std::vector<std::string> slash_spliter(const std::string& str){
    std::vector<std::string> splited;
    std::string element;
    for (auto x : str){
        if (x == '/')
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
void read(const char *file_in, std::map<std::string, std::vector<Point>> &v,
          std::map<std::string, std::string>& constr,
          std::map<std::string, unsigned int>& floors ) {
    //the file needs to be of the following format starting from column 0:
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ column 1 -> id (string)
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ column 8 -> year of construction (string)
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ column 12 -> roof's height (string)
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ column 13 -> ground's height (string)

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
        //auto q = spliter2(splited_line[19]);
        for (unsigned int i = 16; i<splited_line.size(); i=i+2 ) {//Start from 16th element (from wich the coordinate start and continue with step 2

            //Get X coordinate from current line
            size_t j = 0;
            for ( ; j < splited_line[i].length(); j++ ){ if ( isdigit(splited_line[i][j]) ) break; }
            splited_line[i] = splited_line[i].substr(j, splited_line[i].length() - j );
            double id_x = atof(splited_line[i].c_str());

            //Get Y coordinate from next line
            j=0;
            for ( ; j < splited_line[i+1].length(); j++ ){ if ( isdigit(splited_line[i+1][j]) ) break; }
            splited_line[i+1] = splited_line[i+1].substr(j, splited_line[i+1].length() - j );
            double id_y = atof(splited_line[i+1].c_str());

            //Get Z value from ground's height
            double z = stod(splited_line[13]);
            //Get Z value from roof's height
            double roofz = stod(splited_line[12]);


            verts.push_back(Point(id_x,id_y,z,roofz));
        }
        if (line == ""){continue;}
        std::string id = splited_line[1];
        v[id] = verts; //Map vertices of building
        constr[id]=slash_spliter(splited_line[8])[0]; //Map construction year of building
        floors[id]=floor(verts.back().z_r/3); //Map number of floors of building *assume one floor is 3m high*
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

void building_mapper(std::map<std::string, std::map<unsigned int, std::vector<Point>>>& b,std::map<std::string, std::vector<Point>>& v){

    for (auto const& i : v){
        std::map<unsigned int, std::vector<Point>> second;
        std::vector<Point> CW_vert;
        int reverse;
        for (unsigned int j = 0; j < i.second.size()+2; j++) { //faces

            std::vector<Point> vertices;
            if (j == 0) { //First face -> ground face (CW)
                for (int k = i.second.size()-1; k >= 0; k--) {
                    vertices.push_back(i.second[k]);
                }
                for (unsigned int k = 0; k < i.second.size(); k++) {
                    CW_vert.push_back(i.second[k]);
                }
            }
            else if ( j == 1){ //Second face -> roof face (CCW)

                for (unsigned int k = 0; k < i.second.size(); k++) {
                    vertices.push_back(Point(i.second[k].x,i.second[k].y,i.second[k].z_r,i.second[k].z));
                }
            }



            else{
                if(j-2+1==i.second.size()){break;}
                    auto bottom1 = CW_vert[j-2];
                    auto bottom2 = CW_vert[j-2+1];
                    auto up1 = b[i.first].find(1)->second[j-2];
                    auto up2 = b[i.first].find(1)->second[j-2+1];

                    //CCW wall orientation
                vertices.push_back(bottom1);
                vertices.push_back(bottom2);
                vertices.push_back(up2);
                vertices.push_back(up1);


            }
            second[j] = vertices;
            b[i.first] = second;

        }
        b[i.first].find(0)->second.pop_back();
        b[i.first].find(1)->second.pop_back();
    }
}

int main(int argc, const char * argv[]) {
    const char *file_in = "/home/konstantinos/Desktop/TUDelft-Courses/Q3/GEO1004/hw3/data/testCSV.csv";
    const char *file_out = "/home/konstantinos/Desktop/TUDelft-Courses/Q3/GEO1004/hw3/data/test.json";
    std::map<std::string, std::vector<Point>> vertices;
    std::map<std::string, std::map<unsigned int, std::vector<Point>>> buildings;
    std::map<std::string, std::string> const_year;
    std::map<std::string, unsigned int> storeys;
//TODO: Flip ground orientation
//TODO: Subsample (especially with holes) for testing
//TODO: Identify holes in read function where "[[x"

    //Read CSV file, mapping vertices to their building
    read(file_in,vertices,const_year, storeys);

    //Map-> building : { id : building_id , faces : { id :  face_id, vertices : [[v1,...vn]] }}
    building_mapper(buildings,vertices);

    std::map<Point, int> Vertice_mapper;
    std::vector<Point> ordered_verts;
    int count=0;
    for (auto const& b : buildings){
        for (auto const& f : b.second){
            for (auto const& v : f.second){
                if (std::find(ordered_verts.begin(), ordered_verts.end(), v) !=
                    ordered_verts.end()) {
                    continue;}
                ordered_verts.push_back(v);
                Vertice_mapper[v] = count;
                count++;
            }
        }
    }



    //auto b = vertices["0503100000020070"];
    //auto a= buildings["0503100000020070"];



    //Export to CityJson
    std::map<std::string,std::string> cityjson;
    std::fstream fl;
    fl.open(file_out, std::fstream::in | std::fstream::out | std::fstream::trunc);
    fl << "{\n"<<
          "\t\"type\": \"CityJSON\",\n"<<
          "\t\"version\": \"1.0\",\n" <<
          "\t\"CityObjects\": {\n";

    for (auto const& b : buildings){
        fl << "\t\""<<b.first<<"\": {\n";
        fl << "\t\t\"type\": \"Building\",\n";
        fl << "\t\t\"attributes\": {\n";
        fl << "\t\t\t\"yearOfConstruction\": " << const_year[b.first]<<",\n";
        fl << "\t\t\t\"measuredHeight\": " << b.second.at(0).at(0).z_r<<",\n";
        fl << "\t\t\t\"storeysAboveGround\": " << storeys[b.first]<<"},\n";

        fl << "\t\t\"geometry\": [{\n";
        fl << "\t\t\t\"type\": \"Solid\",\n";
        fl << "\t\t\t\"lod\": "<<1.2<<",\n";\
        fl << "\t\t\t\"boundaries\": [[\n";
        for (auto const& f : b.second) {
            fl << "\t\t\t[[";
            for (auto const &v : f.second){
                if (v==f.second.back()){
                    fl << Vertice_mapper[v];
                }
                else fl << Vertice_mapper[v]<<',';
            }
            fl <<"]],\n";
        }
        fl <<"]]}]},\n";
    };
    fl << "},\n\"vertices\": [\n";
    for (auto const& v : ordered_verts) {
        if (v == ordered_verts.back()) {
            fl << "[" << v.x << "," << v.y << "," << v.z << "]\n";
        }else fl << "[" << v.x << "," << v.y << "," << v.z << "],\n";
    }
    fl << "]\n";

    fl <<"}";
    fl.close();
  return 0;
}
