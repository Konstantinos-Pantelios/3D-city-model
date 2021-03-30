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
        unsigned int hole = 0;
        for (unsigned int i = 16; i<splited_line.size(); i=i+2 ) {//Start from 16th element (from wich the coordinate start and continue with step 2

            if(splited_line[i-1].substr( splited_line[i-1].length() - 2 ) == "]]" && splited_line[i].substr(0,2) == "[["){
                verts.pop_back(); //removes the last "closing" vertex as is redundant. -> Polygon(v1,v2,v3,v4,v1)
                hole++; //Set the hole number that is going to be assigned to the rest of the vertices until other hole is found.
            }
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


            verts.push_back(Point(id_x,id_y,z,roofz,hole));
        }
        verts.pop_back(); //removes the last "closing" vertex as is redundant. -> Polygon(v1,v2,v3,v4,v1)
        if (line == ""){continue;}
        std::string id = splited_line[1];
        v[id] = verts; //Map vertices of building
        constr[id]=slash_spliter(splited_line[8])[0]; //Map construction year of building
        floors[id]=floor((verts.back().z_r-verts.back().z)/3); //Map number of floors of building *assume one floor is 3m high*
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
    Point bottom1, bottom2, up1, up2;


    for (auto const& i : v){
        unsigned int curr_hole = 0;
        unsigned int c_hole_start =0;

        std::map<unsigned int, std::vector<Point>> second;
        std::vector<Point> CCW_vert;
        std::vector<Point> temp_roof;
        int reverse;
        for (unsigned int j = 0; j < i.second.size()+2; j++) { //faces -> vertices + 2 (floor&roof) = number of faces

            std::vector<Point> vertices;
            if (j == 0) { //First face -> ground face (CW)
                for (int k = i.second.size()-1; k >= 0; k--) { //CW orientated face
                    //if (i.second[k].hl!=0){continue;} // Checking whether a hole vertex exists and skips it
                    vertices.push_back(i.second[k]);
                }
                for (unsigned int k = 0; k < i.second.size(); k++) { //temporary CCW orientated ground face to set the walls.
                    //if (i.second[k].hl!=0){break;} // Checking whether a hole vertex exists and skips it
                    CCW_vert.push_back(i.second[k]);
                }
            }
            else if ( j == 1){ //Second face -> roof face (CCW)

                for (unsigned int k = 0; k < i.second.size(); k++) {
                    temp_roof.push_back(Point(i.second[k].x,i.second[k].y,i.second[k].z_r,i.second[k].z,i.second[k].hl));
                    //if (i.second[k].hl!=0){ continue;} // Checking whether a hole vertex exists and skips it
                    vertices.push_back(Point(i.second[k].x,i.second[k].y,i.second[k].z_r,i.second[k].z,i.second[k].hl));
                }
            }



            else{
                    if (j - 2 + 1 == i.second.size()) { // Last vertex connect to the front one
                        bottom1 = CCW_vert[j - 2];
                        bottom2 = CCW_vert[c_hole_start];
                        up1 = temp_roof[j - 2];
                        up2 = temp_roof[c_hole_start];
                    } else {
                        if (i.second[j-2+1].hl!=curr_hole){
                            bottom1 = CCW_vert[j - 2];
                            up1 = temp_roof[j - 2];
                            up2 = temp_roof[c_hole_start];
                            bottom2 = CCW_vert[c_hole_start];
                            curr_hole++; //set current hole number
                            c_hole_start=j-2+1;} //set the starting vertex of current hole}
                        else {
                            bottom1 = CCW_vert[j - 2];
                            up1 = temp_roof[j - 2];
                            up2 = temp_roof[j - 2 + 1];
                            bottom2 = CCW_vert[j - 2 + 1];}
                    }

                    //CCW wall orientation
                    vertices.push_back(bottom1);
                    vertices.push_back(bottom2);
                    vertices.push_back(up2);
                    vertices.push_back(up1);
                }


            second[j] = vertices;
            b[i.first] = second;

        }
        //b[i.first].find(0)->second.pop_back();
        //b[i.first].find(1)->second.pop_back();
    }
}

void cjioDS(std::vector<std::vector<Face>>& all,
            std::map<std::string, std::map<unsigned int, std::vector<Point>>>& buildings,
            std::map<Point, int>& Vertice_mapper,
            std::vector<Point> ordered_verts){

    for (auto const& b : buildings){
        Face faces;
        std::vector<Face> build;
        unsigned int at_face_x=0;

        for (auto const& f : b.second) {
            Face faces;
            unsigned int curr_hole=f.second.front().hl;
            std::vector<std::vector<unsigned int>> holes;
            std::vector<unsigned int> temp; //Initialize vector for conditional use at line 290
            bool lever=false;
            for (auto const &v : f.second){
                if (v.hl==curr_hole) {
                    temp.push_back(Vertice_mapper[v]);
                }
                else {
                    holes.push_back(temp);
                    curr_hole=v.hl;
                    temp={};
                    temp.push_back(Vertice_mapper[v]);
                }
            }
            holes.push_back(temp);
            if (at_face_x==0){
                at_face_x++;
                for (int i = holes.size()-1; i >= 0; i--){
                    if (holes.size()>1) {
                        if (i == holes.size()-1) { faces.Exterior = holes[i]; }
                        else { faces.Interior.push_back(holes[i]); }
                    }else {faces.Exterior = holes[i];}
                }
            }
            else {
                for (unsigned int  i = 0; i < holes.size(); i++){
                    if (holes.size()>1) {
                        if (i == 0) { faces.Exterior = holes[i]; }
                        else { faces.Interior.push_back(holes[i]); }
                    }else {faces.Exterior = holes[i];}
                    if (ordered_verts[holes[i].front()].hl!=0 && holes.size()==1){faces.hole_no=ordered_verts[holes[i].front()].hl;}
                }}
            build.push_back(faces);
        }
        all.push_back(build);
    }
}

int main(int argc, const char * argv[]) {
    const char *file_in = "/home/konstantinos/Desktop/TUDelft-Courses/Q3/GEO1004/hw3/data/CSV/tudcampus.csv";
    const char *file_out = "/home/konstantinos/Desktop/TUDelft-Courses/Q3/GEO1004/hw3/data/test.json";
    std::map<std::string, std::vector<Point>> vertices;
    std::map<std::string, std::map<unsigned int, std::vector<Point>>> buildings;
    std::map<std::string, std::string> const_year;
    std::map<std::string, unsigned int> storeys;


//TODO: Adjust "values" to holes


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

    std::vector<std::vector<Face>> all;
    cjioDS(all, buildings, Vertice_mapper, ordered_verts);






    //Export to CityJson
    std::map<std::string,std::string> cityjson;
    std::fstream fl;
    fl.open(file_out, std::fstream::in | std::fstream::out | std::fstream::trunc);
    fl << "{\n"<<
          "\t\"type\": \"CityJSON\",\n"<<
          "\t\"version\": \"1.0\",\n" <<
          "\t\"metadata\": { \"referenceSystem\": \"urn:ogc:def:crs:EPSG::7415\"},\n"
          "\t\"CityObjects\": {\n";

    int b_counter=0;
    for (auto const& b : buildings){
        if(b.first=="0503100000000131"){
            int a = 0;
        }
        int f_counter=0;
        int hole_n=-1;
        unsigned int curr_hole=99999;

        fl << "\t\""<<b.first<<"\": {\n";
        fl << "\t\t\"type\": \"Building\",\n";
        fl << "\t\t\"attributes\": {\n";
        fl << "\t\t\t\"yearOfConstruction\": " << const_year[b.first]<<",\n";
        fl << "\t\t\t\"measuredHeight\": " << b.second.at(0).at(0).z_r -  b.second.at(0).at(0).z<<",\n";
        fl << "\t\t\t\"storeysAboveGround\": " << storeys[b.first]<<"},\n";

        fl << "\t\t\"geometry\": [{\n";
        fl << "\t\t\t\"type\": \"Solid\",\n";
        fl << "\t\t\t\"lod\": "<<1.2<<",\n";\
        fl << "\t\t\t\"boundaries\": [\n"; //Open bracket of Solid

        if (all[b_counter].front().Interior.empty()){ //No-hole case
            fl<<"\t\t\t[\n"; //Open bracket of EXTERIOR SHELL of solid

            for (auto const& f : all[b_counter]) {
                fl<<"\t\t\t\t[["; //MultiSurface brackets - Always double at No-hole case as there wont be any interior boundaries.
                for (auto const& v : f.Exterior){
                    if (v==f.Exterior.back()){fl<<v;} //Last vertex indicex of face does not need a comma
                    else fl<<v<<",";} //Intermediary vertec indices of face DOES need a comma
                if (f.Exterior == all[b_counter].back().Exterior) {
                    fl<<"]]\n"; //MultiSurface brackets - Always double at No-hole case as there wont be any interior boundaries.
                fl<<"\t\t\t\t\t\t\t\t]\n";} //Last face of building needs to close the solids exterior boundary
                else fl<<"]],\n"; ////MultiSurface brackets - Always double at No-hole case as there wont be any interior boundaries.
            }
        }
        else {
            fl << "\t\t\t[\n"; //Open bracket of EXTERIOR SHELL of solid

            auto hole_num = 0;
            bool lever1 = false;
            bool lever2 = false;
            for (auto const &f : all[b_counter]) {



///Writting the Exterior boundary of the Multisurface
  /*              if (all[b_counter][f_counter-1].hole_no != hole_num ){
                    hole_num = all[b_counter][f_counter+1].hole_no;
                    lever1 = true;}
                else lever1 = false;

                if (all[b_counter][f_counter+1].hole_no != hole_num ){
                    hole_num = all[b_counter][f_counter+1].hole_no;
                    lever2 = true;}
                else lever2 = false;

                if (lever1){fl<< "\t\t\t[\n\t\t\t\t[[";}
                else {fl<< "\t\t\t\t[[";}
*/
                fl<<"\t\t\t[["; //Open bracket of EXTERIOR SHELL of solid
                for (auto const &v : f.Exterior) {
                    if (v == f.Exterior.back() && !f.Interior.empty()) { fl << v <<"],"; } //Last vertex indicex of face does not need a comma
                    else if (v == f.Exterior.back() && f.Interior.empty()) {
                        if (f.Exterior == all[b_counter].back().Exterior) { fl << v << "]]\n"; }
                        else if (lever2) { fl << v << "]]\n\t\t\t\t\t\t\t\t\t],\n"; lever2=false; }
                        else { fl << v << "]],\n"; }
                    }
                    else fl << v << ",";} //Intermediary vertec indices of face DOES need a comma

//Writting the Interior boundary of the Multisurface
                for (auto const &in_mulit : f.Interior) {
                    fl<<"[";
                    for (auto const &in : in_mulit) {
                        if (in == in_mulit.back()) {fl << in; } //Last vertex indicex of face does not need a comma
                        else fl << in << ",";} //Intermediary vertec indices of face DOES need a comma
                    if (in_mulit == f.Interior.back()) {fl<< "]],\n";}
                    else fl << "],\n";}


                if (f.Exterior == all[b_counter].back().Exterior) {
                    fl<<"\t\t\t\t\t\t\t\t\t]\n";}

                f_counter++;
            }
        }

        fl<<"\t\t\t],\n";
        fl<<"\t\t\t\"semantics\": {\n";
        fl<<"\t\t\t\t\"surfaces\": [\n";
        fl<<"\t\t\t\t\t{\"type\": \"GroundSurface\"},\n";
        fl<<"\t\t\t\t\t{\"type\": \"WallSurface\"},\n";
        fl<<"\t\t\t\t\t{\"type\": \"RoofSurface\"}],\n";
        fl<<"\t\t\t\t\"values\": [[";
        for (unsigned int c = 0; c < b.second.size(); c++){
            if (c==0){fl<<0<<",";}
            else if (c==1){fl<<2<<",";}
            else {  if (c == b.second.size() - 1) { fl << 1; }
                    else fl << 1 << ",";}
        }

        if (b.first == buildings.rbegin()->first){
            fl <<"]]\n";
            fl<<"\t\t\t\t}\n";
            fl<<"\t\t\t}]\n";
            fl<<"\t\t}\n";}
        else{fl <<"]]\n";
            fl<<"\t\t\t\t}\n";
            fl<<"\t\t\t}]\n";
            fl<<"\t\t},\n";}
        b_counter++;
    };
    fl << "\t},\n\"vertices\": [\n";
    for (auto const& v : ordered_verts) {
        if (v == ordered_verts.back()) { //"fixed" is used to avoid unexpected rounding of the coordinate values
            fl << "[" << v.x<<std::fixed << "," << v.y<<std::fixed << "," << v.z<<std::fixed << "]\n";
        }else fl << "[" << v.x <<std::fixed<< "," << v.y<<std::fixed << "," << v.z<<std::fixed << "],\n";
    }
    fl << "]\n";
    fl <<"}";

    fl.close();
  return 0;
}
