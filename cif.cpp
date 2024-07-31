#include <bits/stdc++.h>
#include <sys/time.h>
#include <sys/stat.h>
#include <dirent.h>

#include <sys/types.h> // required for stat.h

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif
#define SSTR( x ) dynamic_cast< std::ostringstream & >( ( std::ostringstream() << std::dec << x ) ).str()


using namespace std;

double getTime() {
    timeval t;
    gettimeofday(&t,NULL);
    return 1e-6*t.tv_usec + t.tv_sec;
}
inline bool fileExists(const string& name){
	struct stat buffer;
	return (stat (name.c_str(), &buffer) == 0);
}
bool dirExists(const char *path)
{
    struct stat info;
    if(stat( path, &info ) != 0)
        return false;
    else if(info.st_mode & S_IFDIR)
        return true;
    else
        return false;
}
bool createFolder(const string& folderName) {
    if (dirExists(folderName.c_str())) return true;
    mode_t nMode = 0733; // UNIX style permissions
    int nError = 0;
    nError = mkdir(folderName.c_str());
    return nError == 0;
}

string trim(const string& str, const string& whitespace=" \t\r\n") {
    size_t strBegin = str.find_first_not_of(whitespace);
    if (strBegin == string::npos) return "";
    size_t strEnd = str.find_last_not_of(whitespace);
    size_t strRange = strEnd - strBegin + 1;
    return str.substr(strBegin, strRange);
}
vector <string> splitBy(const string& text, char by) {   // split string by by
	vector <string> vys;
	stringstream ss(text);
    string word;
    while(getline(ss, word, by)){
        vys.push_back(word);
    }
    return vys;
}
vector <string> splitByMultiple(const string& text, const string& delims) {
	vector <string> vys;
    stringstream stringStream(text);
    string line;
    while(getline(stringStream, line))
    {
        size_t prev = 0, pos;
        while ((pos = line.find_first_of(delims, prev)) != string::npos) {
            if (pos > prev) vys.push_back(line.substr(prev, pos - prev));
            prev = pos + 1;
        }
        if (prev < line.length()) vys.push_back(line.substr(prev, string::npos));
    }
    return vys;
}
string int2len(int v, int l) {
	string ret = SSTR(v);
	int dig = floor(log10(v + 0.1)) + 1;
	while(ret.length() < l) ret = " " + ret;
	return ret;
}
bool startsWith(string text, string prefix) {
    return text.rfind(prefix, 0) == 0;
}

void writeVector(ostream& out, vector<double>& data, string ending = "\n") {
    for (int i = 0; i < data.size(); i++) {
        if (i != 0) out << ",";
        out << data[i];
    }
    out << ending;
}

template <class T> T ***Create3D(int N1, int N2, int N3)
{
    T *** array = new T ** [N1];
    array[0] = new T * [N1*N2];
    array[0][0] = new T [N1*N2*N3];
    int i,j,k;
    for( i = 0; i < N1; i++) {
        if (i < N1 -1 ) {
            array[0][(i+1)*N2] = &(array[0][0][(i+1)*N3*N2]);
            array[i+1] = &(array[0][(i+1)*N2]);
        }
        for( j = 0; j < N2; j++) {
            if (j > 0) array[i][j] = array[i][j-1] + N3;
        }
    }
    return array;
};

template <class T> void Delete3D(T ***array) {
    delete[] array[0][0];
    delete[] array[0];
    delete[] array;
};

class Double3D {
  public:
    double x;
    double y;
    double z;
    Double3D(){}
    Double3D(double pX, double pY, double pZ) : x(pX), y(pY), z(pZ) {}
    Double3D(vector<string>& parts) {
        x = atof(parts[0].c_str());
        y = atof(parts[1].c_str());
        z = atof(parts[2].c_str());
    }
    Double3D operator + (Double3D const &other) {
         return Double3D(x + other.x, y + other.y, z + other.z);
    }
    Double3D operator - (Double3D const &other) {
         return Double3D(x - other.x, y - other.y, z - other.z);
    }
    Double3D operator * (double const c) {
         return Double3D(c * x, c * y, c * z);
    }
    // dot product
    double operator * (Double3D const &other) {
         return x * other.x + y * other.y + z * other.z;
    }
    // cross product
    Double3D operator ^ (Double3D const &other) {
         return Double3D(y * other.z - z * other.y, z * other.x - x * other.z, x * other.y - y * other.x);
    }
    double distanceTo(Double3D const & other) {
        double dx = x - other.x;
        double dy = y - other.y;
        double dz = z - other.z;
        return sqrt(dx * dx + dy * dy + dz * dz);
    }
};
ostream& operator << (ostream& out, const Double3D& v) {
    out << "(" << v.x << ", " << v.y << ", " << v.z << ")";
    return out;
}
double edgeAngle(Double3D& a, Double3D& b, Double3D& p, Double3D& q) {
    Double3D ab = b - a;
    Double3D ap = p - a;
    Double3D aq = q - a;
    double normABsq = ab * ab;
    double t = (ap * ab) / normABsq;
    Double3D v = ap - ab * t;
    double u = (aq * ab) / normABsq;
    Double3D w = aq - ab * u;

    Double3D normAQB = aq ^ ab;
    Double3D outPoint = q + normAQB;
    double outValue = outPoint * normAQB;
    double onValue = q * normAQB;
    double testedValue = p * normAQB;
    int sign = (outValue - onValue) * (testedValue - onValue) < 0 ? -1 : 1;
    double angle = acos((v * w) / (sqrt(v * v) * sqrt(w * w)));
    if (isnan(angle)) angle = 0;
    return sign * (M_PI - angle);
}
class Triangle {
  public:
    int vertex[3];
};
class Edge {
  public:
    int vertex[2];
    int thirdVertex[2];
};
class ExtraData {
  public:
    double *** & data;
    int length;
    int meanIndexShift;
    int devIndexShift;
    string action;
    string description;
    bool isBasicType;
    ExtraData(double *** & pData, int pLength, string pAction, string pDescription, int pShift, int pDevShift) :
        data(pData), length(pLength), action(pAction), description(pDescription), meanIndexShift(pShift), devIndexShift(pDevShift)
        {
            isBasicType = true;
        }
    ExtraData(double *** & pData, string pDescription) :
        data(pData), description(pDescription)
        {
            isBasicType = false;
        }
};

class Simulation {
  public:
    string shortcut;
    string folder;
    string subfolder;
    int cellCount = 0;
    int cellTypeCount = 0;
    vector<int> cellsPerType;
    vector<double> cellTypeElasticity;
    vector<int> cellTypes;
    vector<string> dataDescriptions;
    unordered_map<string, int> descriptionIndex;
    int stepSize = 0;
    int stepCount = 0;
    int nodes = 0;
    int dataColumnCount = 0;
    bool dataFileContainsData = false;
    double ***data;
    Double3D ***vtkData;
    double ***edgeAngles;
    double ***nodeAngles;
    double ***edgeLengths;
    double ***edgeAnglesDelta;
    double ***nodeAnglesDelta;
    double ***edgeLengthsDelta;
    double ***edgeAnglesDeltaAbs;
    double ***nodeAnglesDeltaAbs;
    double ***edgeLengthsDeltaAbs;
    vector<Triangle> triangles;
    vector<Edge> edges;
    vector<vector<int>> nodeEdges;
    bool dataAllocated = false;
    bool vtkDataAllocated = false;
    bool angleDataAllocated = false;
    vector<ExtraData> extraData;
    vector<int> mainCircleNodes;
    int hollowNodes[2];
    Simulation(string pShortcut, string pFolder, string pSubfolder) : shortcut(pShortcut), folder(pFolder), subfolder(pSubfolder) {
    }
    bool registerExtraData() {
        extraData.push_back(ExtraData(edgeAngles, edges.size(), "mean", "edge_angle_mean", 0, 1));
        extraData.push_back(ExtraData(edgeAngles, edges.size(), "dev", "edge_angle_deviation", -1, 0));
        extraData.push_back(ExtraData(edgeAngles, edges.size(), "skew", "edge_angle_skewness", -2, -1));
        extraData.push_back(ExtraData(nodeAngles, nodes, "mean", "node_angle_mean", 0, 1));
        extraData.push_back(ExtraData(nodeAngles, nodes, "dev", "node_angle_deviation", -1, 0));
        extraData.push_back(ExtraData(nodeAngles, nodes, "skew", "node_angle_skewness", -2, -1));
        extraData.push_back(ExtraData(edgeLengths, edges.size(), "mean", "edge_length_mean", 0, 1));
        extraData.push_back(ExtraData(edgeLengths, edges.size(), "dev", "edge_length_deviation", -1, 0));
        extraData.push_back(ExtraData(edgeLengths, edges.size(), "skew", "edge_length_skewness", -2, -1));
        extraData.push_back(ExtraData(data, "cell_axis_length"));
        extraData.push_back(ExtraData(data, "equator_diameter"));
        extraData.push_back(ExtraData(data, "equator_diameter_min"));
        extraData.push_back(ExtraData(data, "span_x"));
        extraData.push_back(ExtraData(data, "span_y"));
        extraData.push_back(ExtraData(data, "span_z"));
        extraData.push_back(ExtraData(data, "vel_diff_x"));
        extraData.push_back(ExtraData(data, "vel_diff_y"));
        extraData.push_back(ExtraData(data, "vel_diff_z"));
        extraData.push_back(ExtraData(edgeAnglesDelta, edges.size(), "mean", "edge_angle_delta_mean", 0, 1));
        extraData.push_back(ExtraData(edgeAnglesDelta, edges.size(), "dev", "edge_angle_delta_deviation", -1, 0));
        extraData.push_back(ExtraData(edgeAnglesDelta, edges.size(), "skew", "edge_angle_delta_skewness", -2, -1));
        extraData.push_back(ExtraData(nodeAnglesDelta, nodes, "mean", "node_angle_delta_mean", 0, 1));
        extraData.push_back(ExtraData(nodeAnglesDelta, nodes, "dev", "node_angle_delta_deviation", -1, 0));
        extraData.push_back(ExtraData(nodeAnglesDelta, nodes, "skew", "node_angle_delta_skewness", -2, -1));
        extraData.push_back(ExtraData(edgeLengthsDelta, edges.size(), "mean", "edge_length_delta_mean", 0, 1));
        extraData.push_back(ExtraData(edgeLengthsDelta, edges.size(), "dev", "edge_length_delta_deviation", -1, 0));
        extraData.push_back(ExtraData(edgeLengthsDelta, edges.size(), "skew", "edge_length_delta_skewness", -2, -1));
        extraData.push_back(ExtraData(edgeAnglesDeltaAbs, edges.size(), "mean", "edge_angle_delta_abs_mean", 0, 1));
        extraData.push_back(ExtraData(edgeAnglesDeltaAbs, edges.size(), "dev", "edge_angle_delta_abs_deviation", -1, 0));
        extraData.push_back(ExtraData(edgeAnglesDeltaAbs, edges.size(), "skew", "edge_angle_delta_abs_skewness", -2, -1));
        extraData.push_back(ExtraData(nodeAnglesDeltaAbs, nodes, "mean", "node_angle_delta_abs_mean", 0, 1));
        extraData.push_back(ExtraData(nodeAnglesDeltaAbs, nodes, "dev", "node_angle_delta_abs_deviation", -1, 0));
        extraData.push_back(ExtraData(nodeAnglesDeltaAbs, nodes, "skew", "node_angle_delta_abs_skewness", -2, -1));
        extraData.push_back(ExtraData(edgeLengthsDeltaAbs, edges.size(), "mean", "edge_length_delta_abs_mean", 0, 1));
        extraData.push_back(ExtraData(edgeLengthsDeltaAbs, edges.size(), "dev", "edge_length_delta_abs_deviation", -1, 0));
        extraData.push_back(ExtraData(edgeLengthsDeltaAbs, edges.size(), "skew", "edge_length_delta_abs_skewness", -2, -1));
        for(ExtraData& extra: extraData) {
            if(descriptionIndex.count(extra.description) == 0) {
                descriptionIndex[extra.description] = dataDescriptions.size();
                dataDescriptions.push_back(extra.description);
            } else {
                cerr << "    Description '" << extra.description << "' occurs multiple times!" << endl;
                return false;
            }
        }
        return true;
    }
    double calculateAction(double* sourceData, int length, string action, double mean = 0, double dev = 0) {
        double result = 0;
        if (action == "mean") {
            for (int i = 0; i < length; i++) {
                result += sourceData[i];
            }
            return result / length;
        } else if (action == "dev") {
            for (int i = 0; i < length; i++) {
                result += (sourceData[i] - mean) * (sourceData[i] - mean);
            }
            return sqrt(result / length);
        } else if (action == "skew") {
            for (int i = 0; i < length; i++) {
                result += (sourceData[i] - mean) * (sourceData[i] - mean) * (sourceData[i] - mean);
            }
            return dev == 0 ? 0 : result / length / (dev * dev * dev);
        }
        return 0;
    }
    double calculateOneValueOfOneExtraData(ExtraData& extra, int cell, int step, int index) {
        if (extra.isBasicType) return calculateAction(extra.data[cell][step], extra.length, extra.action, data[cell][step][index + extra.meanIndexShift], data[cell][step][index + extra.devIndexShift]);
        // vypocet dalsich dopocitanych features
        if (extra.description == "cell_axis_length") {
           return vtkData[cell][step][hollowNodes[0]].distanceTo(vtkData[cell][step][hollowNodes[1]]);
        }
        if (extra.description == "equator_diameter") {
            double maxDiameter = 0;
            int segmentCount = mainCircleNodes.size() / 2;
            for (int i = 0; i < segmentCount; i++) {
                double diameter = vtkData[cell][step][mainCircleNodes[i]].distanceTo(vtkData[cell][step][mainCircleNodes[i + segmentCount]]);
                if (diameter > maxDiameter) maxDiameter = diameter;
            }
            return maxDiameter;
        }
        if (extra.description == "equator_diameter_min") {
            double minDiameter = 1e100;
            int segmentCount = mainCircleNodes.size() / 2;
            for (int i = 0; i < segmentCount; i++) {
                double diameter = vtkData[cell][step][mainCircleNodes[i]].distanceTo(vtkData[cell][step][mainCircleNodes[i + segmentCount]]);
                if (diameter < minDiameter) minDiameter = diameter;
            }
            return minDiameter;
        }
        for (string dir : {"x", "y", "z"}) {
            if (extra.description == "span_" + dir) {
                return data[cell][step][descriptionIndex[dir + "_max_" + dir]] - data[cell][step][descriptionIndex[dir + "_min_" + dir]];
            }
        }
        for (string dir : {"x", "y", "z"}) {
            if (extra.description == "vel_diff_" + dir) {
                return data[cell][step][descriptionIndex[dir + "_max_vel_" + dir]] - data[cell][step][descriptionIndex[dir + "_min_vel_" + dir]];
            }
        }
        return 0;
    }
    void calculateOneExtraData(ExtraData& extra, int index) {
        for (int i = 0; i < cellCount; i++) {
            for (int j = 0; j < stepCount; j++) {
                data[i][j][index] = calculateOneValueOfOneExtraData(extra, i, j, index);//calculateAction(extra.data[i][j], extra.length, extra.action, data[i][j][index + extra.meanIndexShift], data[i][j][index + extra.devIndexShift]);
            }
        }
    }
    void calculateAllExtraData() {
        for (int i = 0; i < extraData.size(); i++) {
            calculateOneExtraData(extraData[i], dataColumnCount + i);
        }
        //saveData();
        clog << "  Extra data calculated." << endl;
    }
    double angle(int cell, int step, int edge) {
        return edgeAngle(
            vtkData[cell][step][edges[edge].vertex[0]],
            vtkData[cell][step][edges[edge].vertex[1]],
            vtkData[cell][step][edges[edge].thirdVertex[0]],
            vtkData[cell][step][edges[edge].thirdVertex[1]]
        );
    }
    double edgeLength(int cell, int step, int edge) {
        return vtkData[cell][step][edges[edge].vertex[0]].distanceTo(vtkData[cell][step][edges[edge].vertex[1]]);
    }
    string path() {
        return "data/" + folder + "/" + subfolder + "/";
    }
    int dataSize() {
        return cellCount * (dataColumnCount + extraData.size()) * stepCount;
    }
    int vtkDataSize() {
        return cellCount * stepCount * nodes;
    }
    int edgeAngleDataSize() {
        return cellCount * stepCount * edges.size();
    }
    int nodeAngleDataSize() {
        return cellCount * stepCount * nodes;
    }
    int edgeLengthDataSize() {
        return cellCount * stepCount * edges.size();
    }
    string descriptionFile() {
        return path() + "rbc_description_" + subfolder + ".dat";
    }
    string dataFile(int cellNumber) {
        if (dataFileContainsData) return path() + "rbc" + SSTR(cellNumber) + "_data_" + subfolder + ".dat";
        else return path() + "rbc" + SSTR(cellNumber) + "_" + subfolder + ".dat";
    }
    string vtkFile(int cellNumber, int step) {
        return path() + "vtk/rbc" + SSTR(cellNumber) + "_" + SSTR(step * stepSize) + ".vtk";
    }
    bool parseIntValue(string line, string lineStarting, int& property) {
        if (startsWith(line, lineStarting)) {
            vector<string> parts = splitBy(line, ':');
            if (parts.size() == 2) {
                property = atoi(trim(parts[1]).c_str());
                return true;
            }
        }
        return false;
    }
    void setMainCircleAndHollowNodes() {
        if (nodes == 374) {
            vector<int> nodes = {
                1, 20, 21, 22, 23, 24, 25, 26, 27, 2, 28, 29, 30, 31, 32, 33, 34, 35,
                3, 43, 44, 45, 46, 47, 48, 49, 50, 4, 58, 59, 60, 61, 62, 63, 64, 65
            };
            mainCircleNodes = nodes;
            hollowNodes[0] = 0;
            hollowNodes[1] = 5;
        } else if (nodes == 642) {
            vector<int> nodes = {
                609, 41, 630, 161, 640, 7, 468, 123, 469, 34, 450, 118, 445, 6, 371, 97,
                358, 27, 335,  89, 327, 4, 435, 114, 436, 31, 437, 115, 441, 5, 602, 154
            };
            mainCircleNodes = nodes;
            hollowNodes[0] = 12;
            hollowNodes[1] = 17;
        }
        else {
            mainCircleNodes.clear();
            hollowNodes[0] = 0;
            hollowNodes[1] = 0;
            clog << "Number of nodes is different from 374 and 642. Cell equator and cell main axis not known.";
        }
    }
    bool loadInfo() {
        string fileName = path() + "sim_info.dat";
        vector<string> parts;
        clog << "  Loafing simulation info..." << endl;
        if (fileExists(fileName)) {
            ifstream in(fileName);
            string line;
            while (!in.eof()) {
                getline(in, line);
                if (!in.eof()) {
                    line = trim(line);
                    if ((parseIntValue(line, "pocet vsetkych buniek:", cellCount) || parseIntValue(line, "pocet zadanych buniek (spolu):", cellCount)) && cellCount > 0) {
                        cellTypes = vector<int>(cellCount, -1);
                    }
                    if (startsWith(line, "pocet buniek typu")) {
                        parts = splitBy(line, ':');
                        if (parts.size() == 2) {
                            cellsPerType.push_back(atoi(trim(parts[1]).c_str()));
                            cellTypeElasticity.push_back(-1);
                            if (cellsPerType.back() > 0) cellTypeCount = cellsPerType.size();
                            if (atoi(splitBy(parts[0], ' ').back().c_str()) != cellsPerType.size()) {
                                cerr << line << endl << "  Unexpected order of lines in info file!" << endl;
                                return false;
                            }
                        }
                    }
                    parseIntValue(line, "n_steps:", stepSize);
                    parseIntValue(line, "nnode_rbc:", nodes);
                    if (startsWith(line, "bunky typu ")) {
                        parts = splitBy(line, ':');
                        if (parts.size() >= 3 && parts[1] == " ks_rbc") {
                            vector<string> subparts = splitBy(parts[0], ' ');
                            if (subparts.size() == 3) {
                                int type = atoi(subparts[2].c_str()) - 1;
                                if (type >= 0 && type < cellsPerType.size()) {
                                    subparts = splitBy(parts[2], ',');
                                    if (subparts.size() > 0) {
                                        double elasticity = atof(subparts[0].c_str());
                                        if (elasticity >= 0 && elasticity <= 1) {
                                            cellTypeElasticity[type] = elasticity;
                                        }
                                    }
                                }
                            }
                        }
                    }
                    if (startsWith(line, "cell 0 -> ")) {
                        parts = splitBy(line, '|');
                        if (parts.size() >= cellCount) {
                            for (int i = 0; i < cellCount; i++) {
                                string part = trim(parts[i]);
                                if (startsWith(part, "cell " + SSTR(i) + " -> type_")) {
                                    vector<string> subparts = splitBy(part, '_');
                                    if (subparts.size() >= 2) {
                                        int type = atoi(subparts[1].c_str()) - 1;
                                        if (type >= 0) cellTypes[i] = type;
                                    }
                                }
                            }
                        }
                    }
                    if (startsWith(line, SSTR(stepCount * stepSize) + ":")) {
                        stepCount++;
                    }
                }
            }
            in.close();
            if (nodes == 0) {
                cerr << "Info file doesn't contain info on number of nodes!" << endl;
                return false;
            } else {
                setMainCircleAndHollowNodes();
            }
            if (cellCount == 0) {
                cerr << "Info file doesn't contain info on number of cells!" << endl;
                return false;
            }
            if (cellTypeCount == 0) {
                cerr << "Info file doesn't contain info on number of cells in each type!" << endl;
                return false;
            }
            if (stepSize == 0) {
                cerr << "Info file doesn't contain info on step size!" << endl;
                return false;
            }
            for (int i = 0; i < cellTypeCount; i++) {
                if (cellTypeElasticity[i] < 0) {
                    cerr << "Unknown elasticity for some cell types!" << endl;
                    return false;
                }
            }
            for (int i = 0; i < cellCount; i++) {
                if (cellTypes[i] < 0 || cellTypes[i] >= cellTypeCount) {
                    cerr << "Type for some cells not known!" << endl;
                    return false;
                }
            }
            clog << "    Simulation info file found and correct." << endl;
            if (checkBasicDataFiles()) {
                clog << "    Data files for all cells found and correct." << endl;
            } else {
                return false;
            }
            if (loadTriangles()) {
                clog << "    Triangles from VTK file loaded." << endl;
            } else {
                return false;
            }
            if (prepareEdges()) {
                clog << "    Edges from triangles generated." << endl;
            } else {
                return false;
            }
            /*if (checkVtkFiles()) {
                clog << "    VTK files for all cells found.                     " << endl;
            } else {
                return false;
            }*/
            if (loadDataDescription()) {
                clog << "    Data description file parsed successfully." << endl;
            } else {
                return false;
            }
            if (dataDescriptions.size() != dataColumnCount) {
                cerr << "Number of columns in data files does not match number of data descriptions!" << endl;
                return false;
            }
            registerExtraData();
            return true;
        }
        else {
            cerr << "  Info file " << fileName << " not found!" << endl;
            return false;
        }
    }
    bool checkBasicDataFiles() {
        set<int> iterations;
        for (int i = 0; i < cellCount; i++) {
            string fileName = dataFile(i);
            string fileName2;
            if (!fileExists(fileName)) {
                dataFileContainsData = true;
                fileName2 = fileName;
                fileName = dataFile(i);
            }
            if (fileExists(fileName)) {
                ifstream in(fileName.c_str());
                string line;
                int expectedNumber = 0;
                while (!in.eof()) {
                    getline(in, line);
                    if (!in.eof()) {
                        vector<string> parts = splitBy(line, ' ');
                        if (!parts.empty()) {
                            if (dataColumnCount == 0) {
                                dataColumnCount = parts.size();
                            } else {
                                if (dataColumnCount != parts.size()) {
                                    cerr << "    Data file " << fileName << " contains line with unexpected number of items!" << endl;
                                    return false;
                                }
                            }
                            int nextNumber = atoi(parts[0].c_str());
                            if (nextNumber == expectedNumber) {
                                expectedNumber += stepSize;
                            } else {
                                cerr << "    Data file " << fileName << " contains unexpected line starting with '" << parts[0] << "'!" << endl;
                                return false;
                            }
                        }
                    }
                }
                iterations.insert(expectedNumber / stepSize);

            } else {
                cerr << "    Data file " << fileName << " missing!" << endl;
                cerr << "    Data file " << fileName2 << " missing!" << endl;
                return false;
            }
        }
        if (iterations.size() != 1) {
            cerr << "    Number of iterations in data files varies!" << endl;
            return false;
        }
        int iterationCount = *iterations.begin();
        if (iterationCount == stepCount) {
            clog << "    Simulation step count matches data files." << endl;
        } else if (iterationCount == stepCount + 1) {
            stepCount++;
            clog << "    Simulation step count incremented by 1 based on data files. New value: " << stepCount << endl;
        } else {
            cerr << "    Simulation step count (" << stepCount << ") doesn't match data files (" << iterationCount << ")!" << endl;
            return false;
        }
        return true;
    }
    bool checkVtkFiles() {
        for (int i = 0; i < cellCount; i++) {
            //string fileName = path() + "vtk/rbc" + SSTR(i) + ".vtk";
            /*if (!fileExists(fileName)) {
                cerr << "  Data file " << fileName << " missing!" << endl;
                return false;
            }*/
            for (int j = 0; j < stepCount; j++) {
                if (j % 100 == 0) clog << "    Checking VTK file for cell " << i << ", step " << j << "\r";
                string fileName = vtkFile(i, j);
                if (!fileExists(fileName)) {
                    cerr << "  Data file " << fileName << " missing!" << endl;
                    return false;
                }
            }
        }
        return true;
    }
    bool loadDataDescription() {
        string fileName = descriptionFile();
        if (fileExists(fileName)) {
            ifstream in(fileName);
            string line;
            getline(in, line);
            vector<string> parts = splitByMultiple(line, " -[");
            int nextIndex = 1;
            int lastUnlabeledIndex = 0;
            bool readingIndices = false;
            for (string& part: parts) {
                if(part == SSTR(nextIndex)) {
                    if (!readingIndices) {
                        readingIndices = true;
                        lastUnlabeledIndex = nextIndex;
                    }
                    dataDescriptions.push_back("");
                    nextIndex++;
                } else {
                    if (readingIndices) {
                        char extra = lastUnlabeledIndex + 1 == nextIndex ? ' ' : 'x';
                        while (lastUnlabeledIndex < nextIndex) {
                            if (lastUnlabeledIndex < 1 || lastUnlabeledIndex > dataDescriptions.size()) {
                                cerr << "  Description file parsing failed!" << endl;
                                return false;
                            }
                            dataDescriptions[lastUnlabeledIndex - 1] = part;
                            if (extra != ' ') {
                                dataDescriptions[lastUnlabeledIndex - 1] += "_";
                                dataDescriptions[lastUnlabeledIndex - 1] += extra;
                                extra++;
                            }
                            lastUnlabeledIndex++;
                        }
                        readingIndices = false;
                    }
                }
            }
            if (!dataDescriptions.empty() && dataDescriptions.back() == "") {
                cerr << "  Not all data descriptions parsed properly!" << endl;
                return false;
            }
            for (int i = 0; i < dataDescriptions.size(); i++) {
                descriptionIndex[dataDescriptions[i]] = i;
            }
        } else {
            cerr << "  Description file " << fileName << " not found!" << endl;
            return false;
        }
        return true;
    }
    void allocateData() {
        data = Create3D<double>(cellCount, stepCount, dataColumnCount + extraData.size());
        dataAllocated = true;
    }
    void deallocateData() {
        Delete3D(data);
        dataAllocated = false;
    }
    void allocateVtkData() {
        vtkData = Create3D<Double3D>(cellCount, stepCount, nodes);
        vtkDataAllocated = true;
    }
    void deallocateVtkData() {
        Delete3D(vtkData);
        vtkDataAllocated = false;
    }
    void allocateAngleData() {
        edgeAngles = Create3D<double>(cellCount, stepCount, edges.size());
        nodeAngles = Create3D<double>(cellCount, stepCount, nodes);
        edgeLengths = Create3D<double>(cellCount, stepCount, edges.size());
        edgeAnglesDelta = Create3D<double>(cellCount, stepCount, edges.size());
        nodeAnglesDelta = Create3D<double>(cellCount, stepCount, nodes);
        edgeLengthsDelta = Create3D<double>(cellCount, stepCount, edges.size());
        edgeAnglesDeltaAbs = Create3D<double>(cellCount, stepCount, edges.size());
        nodeAnglesDeltaAbs = Create3D<double>(cellCount, stepCount, nodes);
        edgeLengthsDeltaAbs = Create3D<double>(cellCount, stepCount, edges.size());
        angleDataAllocated = true;
    }
    void deallocateAngleData() {
        Delete3D(edgeAngles);
        Delete3D(nodeAngles);
        Delete3D(edgeLengths);
        Delete3D(edgeAnglesDelta);
        Delete3D(nodeAnglesDelta);
        Delete3D(edgeLengthsDelta);
        Delete3D(edgeAnglesDeltaAbs);
        Delete3D(nodeAnglesDeltaAbs);
        Delete3D(edgeLengthsDeltaAbs);
        angleDataAllocated = false;
    }
    bool loadTriangles() {
        string fileName = vtkFile(0, 0);
        if (fileExists(fileName)) {
            ifstream in(fileName.c_str());
            string line;
            vector<string> parts;
            int nextIndex = -1;
            while (!in.eof()) {
                getline(in, line);
                if (!in.eof()) {
                    if (startsWith(line, "TRIANGLE_STRIPS")) {
                        parts = splitBy(line, ' ');
                        if (parts.size() >= 2) {
                            int triangleCount = atoi(parts[1].c_str());
                            if (triangleCount > 0) {
                                triangles = vector<Triangle>(triangleCount);
                                nextIndex = 0;
                            }
                        }
                    } else if (nextIndex >= 0 && nextIndex < triangles.size()) {
                        parts = splitBy(line, ' ');
                        if (parts.size() == 4 && parts[0] == "3") {
                            int vertex[3];
                            for (int j = 0; j < 3; j++) {
                                vertex[j] = atoi(parts[1 + j].c_str());
                                if (vertex[j] < 0 || vertex[j] >= nodes) {
                                    cerr << "  Incorrect vertex index in VTK file!" << endl;
                                    return false;
                                }
                                triangles[nextIndex].vertex[j] = vertex[j];
                            }
                            nextIndex++;
                        }
                    }
                }
            }
            if (nextIndex != triangles.size()) {
                cerr << "    Incorrect number of vertex lines in VTK file!" << endl;
                return false;
            }
        } else {
            cerr << "    VTK file " << fileName << " missing!" << endl;
            return false;
        }
        return true;
    }
    int edgeToId(int vertex1, int vertex2) {
        //return min(vertex1, vertex2) * nodes + max(vertex1, vertex2);
        return vertex1 * nodes + vertex2;
    }
    bool prepareEdges() {
        unordered_map<int, int> unpairedEdges;
        nodeEdges = vector<vector<int>>(nodes);
        for (Triangle& triangle: triangles) {
            for (int i = 0; i < 3; i++) {
                int edgeId = edgeToId(triangle.vertex[i], triangle.vertex[(i + 1) % 3]);
                int edgeIdInv = edgeToId(triangle.vertex[(i + 1) % 3], triangle.vertex[i]);
                if (unpairedEdges.count(edgeIdInv) > 0) {
                    edges.push_back(Edge());
                    Edge& edge = edges.back();
                    edge.vertex[0] = triangle.vertex[i];
                    edge.vertex[1] = triangle.vertex[(i + 1) % 3];
                    edge.thirdVertex[0] = triangle.vertex[(i + 2) % 3];
                    edge.thirdVertex[1] = unpairedEdges[edgeIdInv];
                    unpairedEdges.erase(edgeIdInv);
                    nodeEdges[edge.vertex[0]].push_back(edges.size() - 1);
                    nodeEdges[edge.vertex[1]].push_back(edges.size() - 1);
                } else {
                    unpairedEdges[edgeId] = triangle.vertex[(i + 2) % 3];
                }
            }
        }
        if (!unpairedEdges.empty()) {
            cerr << "    Some triangles have unpaired edges!" << endl;
            return false;
        }
        if (edges.size() + 2 != nodes + triangles.size()) {
            cerr << "    Unexpected number of edges (E + 2 != V + F)!" << endl;
            return false;
        }
        return true;
    }
    bool loadDataFromOrigin() {
        if (!dataAllocated) allocateData();
        clog << "  Loading data..." << endl;
        for (int i = 0; i < cellCount; i++) {
            clog << "    Cell " << i << "\r";
            string fileName = dataFile(i);
            if (fileExists(fileName)) {
                ifstream in(fileName.c_str());
                string line;
                int step = 0;
                while (!in.eof()) {
                    getline(in, line);
                    if (!in.eof()) {
                        vector<string> parts = splitBy(line, ' ');
                        if (parts.size() == dataColumnCount) {
                            for (int j = 0; j < dataColumnCount; j++) {
                                data[i][step][j] = atof(parts[j].c_str());
                            }
                            step++;
                        }
                    }
                }
                if (step != stepCount) {
                    cerr << "    Data file " << fileName << " contains unexpected number of lines!" << endl;
                    return false;
                }
            } else {
                cerr << "    Data file " << fileName << " missing!" << endl;
                return false;
            }
        }
        clog << "  Data loaded successfully." << endl;
        return true;
    }
    bool loadDataFromBinary() {
        if (!dataAllocated) allocateData();
        string fileName = path() + "bindata/data.bin";
        if (fileExists(fileName.c_str())) {
            ifstream in(fileName.c_str(), ios::binary);
            for (int i = 0; i < cellCount; i++) for (int j = 0; j < stepCount; j++) {
                in.read((char *)&data[i][j][0], dataColumnCount * sizeof(double));
            }
            in.close();
        } else {
            return false;
        }
        return true;
    }
    bool loadData() {
        if (!loadDataFromBinary()) {
            clog << "  Binary data not found. Loading from origin." << endl;
            if (!loadDataFromOrigin()) return false;
            else saveData();
        } else {
            clog << "  Data loaded successfully from binary file." << endl;
        }
        return true;
    }
    bool saveData() {
        if (createFolder(path() + "bindata")) {
            string fileName = path() + "bindata/data.bin";
           	ofstream out(fileName.c_str(), ios::binary);
            for (int i = 0; i < cellCount; i++) for (int j = 0; j < stepCount; j++) {
                out.write((const char *)&data[i][j][0], dataColumnCount * sizeof(double));
            }
            out.close();
        } else {
            cerr << "  Folder 'bindata/' not created!" << endl;
            return false;
        }
        return true;
    }
    bool loadVtkDataFromOrigin() {
        if (!vtkDataAllocated) allocateVtkData();
        clog << "  Loading VTK data..." << endl;
        vector <string> parts;
        for (int i = 0; i < cellCount; i++) {
            for (int j = 0; j < stepCount; j++) {
                if (j%100 == 0) clog << "    Cell " << i << ", step " << j << "      \r";
                string fileName = vtkFile(i, j);
                if (fileExists(fileName)) {
                    ifstream in(fileName.c_str());
                    string line;
                    int node = -1;
                    while (!in.eof() && node < nodes) {
                        getline(in, line);
                        if (!in.eof()) {
                            if (node < 0 && startsWith(line, "POINTS ")) {
                                parts = splitBy(line, ' ');
                                if (parts.size() > 1) {
                                    int nodeCount = atoi(parts[1].c_str());
                                    if (nodeCount != nodes) {
                                        cerr << "    Incorrect number of nodes in VTK file " << fileName << "!" << endl;
                                        return false;
                                    }
                                    node = 0;
                                }
                            } else if (node >= 0) {
                                vector<string> parts = splitBy(line, ' ');
                                if (parts.size() == 3) {
                                    vtkData[i][j][node] = Double3D(parts);
                                    node++;
                                }
                            }
                        }
                    }
                    if (node != nodes) {
                        cerr << "    VTK file " << fileName << " contains unexpected number of nodes!" << endl;
                        return false;
                    }
                } else {
                    cerr << "    VTK file " << fileName << " missing!" << endl;
                    return false;
                }
            }
        }
        clog << "  VTK data loaded successfully." << endl;
        return true;
    }
    bool loadVtkDataFromBinary() {
        if (!vtkDataAllocated) allocateVtkData();
        string fileName = path() + "bindata/vtk_data.bin";
        if (fileExists(fileName.c_str())) {
            ifstream in(fileName.c_str(), ios::binary);
            in.read((char *)&vtkData[0][0][0], vtkDataSize() * sizeof(Double3D));
            in.close();
        } else {
            return false;
        }
        return true;
    }
    bool loadVtkData() {
        if (!loadVtkDataFromBinary()) {
            clog << "  Binary VTK data not found. Loading from origin." << endl;
            if (!loadVtkDataFromOrigin()) return false;
            else saveVtkData();
        } else {
            clog << "  VTK data loaded successfully from binary file." << endl;
        }
        return true;
    }
    bool saveVtkData() {
        if (createFolder(path() + "bindata")) {
            string fileName = path() + "bindata/vtk_data.bin";
           	ofstream out(fileName.c_str(), ios::binary);
            out.write((const char *)&vtkData[0][0][0], vtkDataSize() * sizeof(Double3D));
            out.close();
        } else {
            cerr << "  Folder 'bindata/' not created!" << endl;
            return false;
        }
        return true;
    }
    bool loadEdgeAngleDataFromOrigin() {
        if (!angleDataAllocated) allocateAngleData();
        clog << "  Generating edge angle data..." << endl;
        for (int i = 0; i < cellCount; i++) {
            for (int j = 0; j < stepCount; j++) {
                if (j%100 == 0) clog << "    Cell " << i << ", step " << j << "      \r";
                for (int k = 0; k < edges.size(); k++) {
                    edgeAngles[i][j][k] = angle(i, j, k);
                    edgeAnglesDelta[i][j][k] = edgeAngles[i][j][k] - edgeAngles[i][0][k];
                    edgeAnglesDeltaAbs[i][j][k] = fabs(edgeAnglesDelta[i][j][k]);
                }
            }
        }
        clog << "  Edge angle data generated successfully." << endl;
        return true;
    }
    bool loadNodeAngleDataFromOrigin() {
        if (!angleDataAllocated) allocateAngleData();
        clog << "  Generating node angle data..." << endl;
        for (int i = 0; i < cellCount; i++) {
            for (int j = 0; j < stepCount; j++) {
                if (j%100 == 0) clog << "    Cell " << i << ", step " << j << "      \r";
                for (int k = 0; k < nodes; k++) {
                    nodeAngles[i][j][k] = 2 * M_PI;
                    for (int e: nodeEdges[k]) {
                        nodeAngles[i][j][k] -= edgeAngles[i][j][e];
                    }
                    nodeAnglesDelta[i][j][k] = nodeAngles[i][j][k] - nodeAngles[i][0][k];
                    nodeAnglesDeltaAbs[i][j][k] = fabs(nodeAnglesDelta[i][j][k]);
                }
            }
        }
        clog << "  Node angle data generated successfully." << endl;
        return true;
    }
    bool loadEdgeLengthDataFromOrigin() {
        if (!angleDataAllocated) allocateAngleData();
        clog << "  Generating edge length data..." << endl;
        for (int i = 0; i < cellCount; i++) {
            for (int j = 0; j < stepCount; j++) {
                if (j%100 == 0) clog << "    Cell " << i << ", step " << j << "      \r";
                for (int k = 0; k < edges.size(); k++) {
                    edgeLengths[i][j][k] = edgeLength(i, j, k);
                    edgeLengthsDelta[i][j][k] = edgeLengths[i][j][k] - edgeLengths[i][0][k];
                    edgeLengthsDeltaAbs[i][j][k] = fabs(edgeLengthsDelta[i][j][k]);
                }
            }
        }
        clog << "  Edge length data generated successfully." << endl;
        return true;
    }
    bool loadEdgeAngleDataFromBinary() {
        if (!angleDataAllocated) allocateAngleData();
        string fileName = path() + "bindata/edge_angle_data.bin";
        if (fileExists(fileName.c_str())) {
            ifstream in(fileName.c_str(), ios::binary);
            in.read((char *)&edgeAngles[0][0][0], edgeAngleDataSize() * sizeof(double));
            in.read((char *)&edgeAnglesDelta[0][0][0], edgeAngleDataSize() * sizeof(double));
            in.read((char *)&edgeAnglesDeltaAbs[0][0][0], edgeAngleDataSize() * sizeof(double));
            in.close();
        } else {
            return false;
        }
        return true;
    }
    bool loadNodeAngleDataFromBinary() {
        if (!angleDataAllocated) allocateAngleData();
        string fileName = path() + "bindata/node_angle_data.bin";
        if (fileExists(fileName.c_str())) {
            ifstream in(fileName.c_str(), ios::binary);
            in.read((char *)&nodeAngles[0][0][0], nodeAngleDataSize() * sizeof(double));
            in.read((char *)&nodeAnglesDelta[0][0][0], nodeAngleDataSize() * sizeof(double));
            in.read((char *)&nodeAnglesDeltaAbs[0][0][0], nodeAngleDataSize() * sizeof(double));
            in.close();
        } else {
            return false;
        }
        return true;
    }
    bool loadEdgeLengthDataFromBinary() {
        if (!angleDataAllocated) allocateAngleData();
        string fileName = path() + "bindata/edge_length_data.bin";
        if (fileExists(fileName.c_str())) {
            ifstream in(fileName.c_str(), ios::binary);
            in.read((char *)&edgeLengths[0][0][0], edgeLengthDataSize() * sizeof(double));
            in.read((char *)&edgeLengthsDelta[0][0][0], edgeLengthDataSize() * sizeof(double));
            in.read((char *)&edgeLengthsDeltaAbs[0][0][0], edgeLengthDataSize() * sizeof(double));
            in.close();
        } else {
            return false;
        }
        return true;
    }
    bool loadAngleData() {
        if (!loadEdgeAngleDataFromBinary()) {
            clog << "  Binary edge angle data not found. Generating from origin." << endl;
            if (!loadEdgeAngleDataFromOrigin()) return false;
            else saveEdgeAngleData();
        } else {
            clog << "  Edge angle data loaded successfully from binary file." << endl;
        }
        if (!loadNodeAngleDataFromBinary()) {
            clog << "  Binary node angle data not found. Generating from origin." << endl;
            if (!loadNodeAngleDataFromOrigin()) return false;
            else saveNodeAngleData();
        } else {
            clog << "  Node angle data loaded successfully from binary file." << endl;
        }
        if (!loadEdgeLengthDataFromBinary()) {
            clog << "  Binary edge length data not found. Generating from origin." << endl;
            if (!loadEdgeLengthDataFromOrigin()) return false;
            else saveEdgeLengthData();
        } else {
            clog << "  Edge length data loaded successfully from binary file." << endl;
        }
        return true;
    }
    bool saveEdgeAngleData() {
        if (createFolder(path() + "bindata")) {
            string fileName = path() + "bindata/edge_angle_data.bin";
           	ofstream out(fileName.c_str(), ios::binary);
            out.write((const char *)&edgeAngles[0][0][0], edgeAngleDataSize() * sizeof(double));
            out.write((const char *)&edgeAnglesDelta[0][0][0], edgeAngleDataSize() * sizeof(double));
            out.write((const char *)&edgeAnglesDeltaAbs[0][0][0], edgeAngleDataSize() * sizeof(double));
            out.close();
        } else {
            cerr << "  Folder 'bindata/' not created!" << endl;
            return false;
        }
        return true;
    }
    bool saveNodeAngleData() {
        if (createFolder(path() + "bindata")) {
            string fileName = path() + "bindata/node_angle_data.bin";
           	ofstream out(fileName.c_str(), ios::binary);
            out.write((const char *)&nodeAngles[0][0][0], nodeAngleDataSize() * sizeof(double));
            out.write((const char *)&nodeAnglesDelta[0][0][0], nodeAngleDataSize() * sizeof(double));
            out.write((const char *)&nodeAnglesDeltaAbs[0][0][0], nodeAngleDataSize() * sizeof(double));
            out.close();
        } else {
            cerr << "  Folder 'bindata/' not created!" << endl;
            return false;
        }
        return true;
    }
    bool saveEdgeLengthData() {
        if (createFolder(path() + "bindata")) {
            string fileName = path() + "bindata/edge_length_data.bin";
           	ofstream out(fileName.c_str(), ios::binary);
            out.write((const char *)&edgeLengths[0][0][0], edgeLengthDataSize() * sizeof(double));
            out.write((const char *)&edgeLengthsDelta[0][0][0], edgeLengthDataSize() * sizeof(double));
            out.write((const char *)&edgeLengthsDeltaAbs[0][0][0], edgeLengthDataSize() * sizeof(double));
            out.close();
        } else {
            cerr << "  Folder 'bindata/' not created!" << endl;
            return false;
        }
        return true;
    }
    void loadAllData() {
        loadData();
        loadVtkData();
        loadAngleData();
    }
    void showInfo(ostream& out) {
        out << "  ---------------" << endl;
        out << "  Simulation info" << endl;
        out << "  ---------------" << endl;
        out << "           Number of cells: " << cellCount << endl;
        out << "      Number of cell types: " << cellTypeCount << endl;
        out << "  Number of cells per type:";
        for (int i = 0; i < cellTypeCount; i++) {
            out << " " << cellsPerType[i];
        }
        out << endl;
        out << "       Elasticity per type:";
        for (int i = 0; i < cellTypeCount; i++) {
            out << " " << cellTypeElasticity[i];
        }
        out << endl;
        out << "             Type per cell:";
        for (int i = 0; i < cellCount; i++) {
            out << " " << cellTypes[i];
        }
        out << endl;
        out << "            Nodes per cell: " << nodes << endl;
        out << "                 Triangles: " << triangles.size() << endl;
        out << "                     Edges: " << edges.size() << endl;
        out << "                 Step size: " << stepSize << endl;
        out << "                Step count: " << stepCount << endl;
        out << "              Data columns: " << dataColumnCount << endl;
    }
    void showDataDescription(ostream& out) {
        out << "   Data column description: " << endl;
        for (int i = 0; i < dataDescriptions.size(); i++) {
            out << "             " << int2len(i, 2) << ": " << dataDescriptions[i] << endl;
        }
    }
    void showDataSizeInfo(ostream& out) {
        out << "                Data size: " << dataSize() << " values" << endl;
        out << "                           " << (dataSize() * 8)/1000000 << " MB (double)" << endl;
    }
    void showVtkDataSizeInfo(ostream& out) {
        out << "            VTK data size: " << (vtkDataSize() * 3) << " values" << endl;
        out << "                           " << (vtkDataSize() * sizeof(Double3D))/1000000 << " MB (double)" << endl;
    }
    void showAngleDataSizeInfo(ostream& out) {
        int valueCount = edgeAngleDataSize() + nodeAngleDataSize() + edgeLengthDataSize();
        out << "          angle data size: " << valueCount << " values" << endl;
        out << "                           " << (valueCount * 8)/1000000 << " MB (double)" << endl;
    }
    bool exportEdgeAngles(int step) {
        if (createFolder(path() + "export")) {
            clog << "  Exporting edge angles..." << endl;
            for (int i = 0; i < cellCount; i++) {
                for (int j = 0; j < stepCount; j += step) {
                    clog << "    cell " << i << " / step " << j << "      \r";
                    string fileName = path() + "export/rbc" + SSTR(i) + "_typ" + SSTR(cellTypes[i]) + "_" + SSTR(j * stepSize) + ".csv";
                    ofstream out(fileName.c_str());
                    out.precision(16);
                    for (int k = 0; k < edges.size(); k++) {
                        out << edgeAngles[i][j][k] << endl;
                    }
                    out.close();
                }
            }
            clog << "  Edge angles exported.                " << endl;
        } else {
            cerr << "  Folder 'export/' not created!" << endl;
            return false;
        }
        return true;
    }
    bool exportSurfaceAndVolume() {
        if (createFolder(path() + "export")) {
            clog << "  Exporting surface and volume info..." << endl;
            int surfaceIndex = descriptionIndex["surface"];
            int volumeIndex = descriptionIndex["volume"];
            if (surfaceIndex == 0) {
                cerr << "  Surface info not present in data!" << endl;
                return false;
            }
            if (volumeIndex == 0) {
                cerr << "  Volume info not present in data!" << endl;
                return false;
            }
            for (int i = 0; i < cellCount; i++) {
                string fileName = path() + "export/surface_and_volume_rbc" + SSTR(i) + "_typ" + SSTR(cellTypes[i]) + ".csv";
                ofstream out(fileName.c_str());
                for (int j = 0; j < stepCount; j++) {
                    if (j%100 == 0) clog << "    cell " << i << " / step " << j << "      \r";
                    out.precision(16);
                    out << j * stepSize << "," << data[i][j][surfaceIndex] << "," << data[i][j][volumeIndex] << endl;
                }
                out.close();
            }
            clog << "  Surface and volume exported.                " << endl;
        } else {
            cerr << "  Folder 'export/' not created!" << endl;
            return false;
        }
        return true;
    }
    bool outputSelectedFeatures(string featureList, bool outputHeader = true, int frames = 0) {
        vector<string> features;
        if (featureList == "all") {
            featureList = "rbc_center_position_x,rbc_center_position_y,rbc_center_position_z,rbc_velocity_x,rbc_velocity_y,rbc_velocity_z,volume,surface";
            features = splitBy(featureList, ',');
            for (auto& extra: extraData) {
                features.push_back(extra.description);
            }
        } else {
            features = splitBy(featureList, ',');
        }
        vector<int> featureIndices;
        int maxIndex = 0;
        for(auto& feature: features) {
            if (descriptionIndex.count(feature) == 0) {
                cerr << "Feature '" << feature << "' not present in data. Output not created." << endl;
                return false;
            } else {
                featureIndices.push_back(descriptionIndex[feature]);
                maxIndex = max(maxIndex, descriptionIndex[feature]);
            }
        }
        clog << "  Outputting selected features..." << endl;
        cout.precision(16);
        if (outputHeader) {
            cout << "cell,step,cell_type";
            for(auto& feature: features) cout << "," << feature;
            if (frames > 1) {
                for(auto& feature: features) cout << "," << feature << "_std";
            }
            cout << endl;
        }
        if (frames <= 1) {
            for (int i = 0; i < cellCount; i++) {
                for (int j = 0; j < stepCount; j++) {
                    if (j % 100 == 0) clog << i << " " << j << "      \r";
                    cout << i << "," << j << "," << cellTypes[i];
                    for(int index: featureIndices) {
                       cout << "," << data[i][j][index];
                    }
                    cout << endl;
                }
            }
        } else {
            for (int i = 0; i < cellCount; i++) {
                for (int j = 0; j < stepCount - frames + 1; j++) {
                    if (j % 100 == 0) clog << i << " " << j << "      \r";
                    cout << i << "," << j << "," << cellTypes[i];
                    vector<double> means(maxIndex + 1, 0);
                    for(int index: featureIndices) {
                        double mean = 0;
                        for (int k = 0; k < frames; k++) {
                            mean += data[i][j + k][index];
                        }
                        mean /= frames;
                        cout << "," << mean;
                        means[index] = mean;
                    }
                    for(int index: featureIndices) {
                        double stdev = 0;
                        double mean = means[index];
                        for (int k = 0; k < frames; k++) {
                            stdev += (data[i][j + k][index] - mean) * (data[i][j + k][index] - mean);
                        }
                        stdev = sqrt(stdev / frames);
                        cout << "," << stdev;
                    }
                    cout << endl;
                }
            }
        }
        clog << "  Selected features output.                " << endl;
        return true;
    }
    bool showFeatureStats(string feature, int frames) {
        if (feature == "x") feature = "rbc_center_position_x";
        if (feature == "y") feature = "rbc_center_position_y";
        if (feature == "z") feature = "rbc_center_position_z";
        int featureIndex = 0;
        if (descriptionIndex.count(feature) == 0) {
            cerr << "Feature '" << feature << "' not present in data. Stats not created." << endl;
            return false;
        } else {
            featureIndex = descriptionIndex[feature];
        }

        double minDelta = 1e30;
        double maxDelta = -1e30;
        double meanDelta = 0;
        int deltaCount = 0;
        for (int i = 0; i < cellCount; i++) {
            for (int j = 0; j < stepCount - frames + 1; j++) {
                double minValue = 1e30;
                double maxValue = -1e30;
                for (int k = 0; k < frames; k++) {
                    minValue = min(minValue, data[i][j + k][featureIndex]);
                    maxValue = max(minValue, data[i][j + k][featureIndex]);
                }
                double delta = maxValue - minValue;
                minDelta = min(minDelta, delta);
                maxDelta = max(maxDelta, delta);
                meanDelta += delta;
                deltaCount++;
            }
        }
        meanDelta /= deltaCount;
        cout << "=======================" << endl;
        cout << "     Feature: " << feature << endl;
        cout << "      Frames: " << frames << endl;
        cout << "   Delta min: " << minDelta << endl;
        cout << "   Delta max: " << maxDelta << endl;
        cout << "  Delta mean: " << meanDelta << endl;
        return true;
    }
};

vector<Simulation> simulations = {
    Simulation("sym1", "sim_from_symmetry2", "sim_9xKSin1_Aa_hct10_seed02_1"),
    Simulation("sym2", "sim_from_symmetry2", "sim_9xKSin1_Aa_rbc54_seed01_1"),
    Simulation("sym3", "sim_from_symmetry2", "sim_9xKSin1_Aa_rbc54_seed02_1"),
    Simulation("C1", "sim_C_stvorcovy_prazdny", "sim_4xKSin1_Ca_rbc36_seed03_1"),
    Simulation("C2", "sim_C_stvorcovy_prazdny", "sim_4xKSin1_Ca_rbc36_seed04_1"),
    Simulation("B1", "sim_B_plytky_prazdny", "sim_Kika_channel_bez_prekazok_4xKs_node374_rbc36_seed01"),
    Simulation("B2", "sim_B_plytky_prazdny", "sim_Kika_channel_bez_prekazok_4xKs_node374_rbc36_seed02"),
    Simulation("B3", "sim_B_plytky_prazdny", "sim_Kika_channel_bez_prekazok_4xKs_node642_rbc36_seed01"),
    Simulation("D1", "sim_D_obdlznikovy_s_prekazkami", "sim_4xKSin1_cyl5_Ca_rbc36_seed05_1"),
    Simulation("D2", "sim_D_obdlznikovy_s_prekazkami", "sim_4xKSin1_cyl5_Ca_rbc36_seed06_1")
};

map<string,string> parseArguments(int argc, char* argv[]) {
    map <string,string> option;
    // default values
    option["export"] = "";
    option["show"] = "";
    option["features"] = "";

    for(int i = 1; i < argc - 1; i++) if (argv[i][0] == '-'){
	    string opt = argv[i];
        while (opt.length() > 0 && opt[0] == '-') opt = opt.substr(1);
        option[opt] = argv[i + 1];
	}
	return option;
}
int parsePositiveInt(string text) {
    int value = atoi(text.c_str());
    if (value > 0) return value;
    return 0;
}

int main(int argc, char* argv[]) {

	double start = getTime();
	string shortcut = argc > 1 ? argv[1] : "";
	if (shortcut == "") {
        cout << "Enter simulation shortcut: ";
        cin >> shortcut;
	}
	map<string,string> option = parseArguments(argc, argv);

	int simId = -1;
	for (int i = 0; i < simulations.size(); i++) {
        if (simulations[i].shortcut == shortcut) {
            simId = i;
            break;
        }
	}
	if (simId < 0) {
        clog << "Simulation with the given shortcut '" << shortcut << "' doesn't exist!";
        return 0;
	}
	Simulation& simulation = simulations[simId];
	if(!simulation.loadInfo()) {
        clog << "Problem during simulation loading occurs!" << endl;
        return 0;
	}
	clog << endl;
	//simulation.showInfo(clog);
	//simulation.showDataDescription(clog);
    //simulation.showDataSizeInfo(clog);
    //simulation.showVtkDataSizeInfo(clog);
    //simulation.showAngleDataSizeInfo(clog);
    simulation.loadAllData();


    simulation.calculateAllExtraData();


    if (option["export"] == "edge_angles") {
        int step = parsePositiveInt(option["step"]);
        if (step > 0) {
            simulation.exportEdgeAngles(step);
        }
    }

    if (option["export"] == "surface_and_volume") {
        simulation.exportSurfaceAndVolume();
    }

    if (option["features"] != "") {
        simulation.outputSelectedFeatures(option["features"], option["header"] != "", parsePositiveInt(option["frames"]));
    }

    if (option["stats"] != "") {
        simulation.showFeatureStats(option["stats"], parsePositiveInt(option["frames"]));
    }

    if (option["info"] == "fields") {
        clog << endl << "  Info:" << endl;
        simulation.showDataDescription(clog);
    }


    /*for (int k = 0; k < simulation.edges.size(); k++) {
        Double3D& a = simulation.vtkData[0][0][simulation.edges[k].vertex[0]];
        Double3D& b = simulation.vtkData[0][0][simulation.edges[k].vertex[1]];
        cout << a.x << "," << a.y << "," << a.z << "," << b.x << "," << b.y << "," << b.z << "," << simulation.edgeAngles[0][0][k] << endl;
    }*/
    /*for (int k = 0; k < simulation.nodes; k++) {
        Double3D& a = simulation.vtkData[0][0][k];
        cout << a.x << "," << a.y << "," << a.z << "," << simulation.nodeAngles[0][0][k] << endl;
    }*/


    /*int i1 = simulation.descriptionIndex["node_angle_mean"];
    int i2 = simulation.descriptionIndex["edge_angle_mean"];
    int i3 = simulation.descriptionIndex["node_angle_deviation"];
    int i4 = simulation.descriptionIndex["edge_angle_deviation"];
    int i5 = simulation.descriptionIndex["node_angle_skewness"];
    int i6 = simulation.descriptionIndex["edge_angle_skewness"];
    */

    /*for (int i = 0; i < simulation.cellCount; i++) {
        for (int j = 0; j < simulation.stepCount; j++) {
            if (j % 100 == 0) clog << i << " " << j << "      \r";
            cout << simulation.cellTypes[i] << ","
                << simulation.data[i][j][i1] << "," << simulation.data[i][j][i2] << ","
                << simulation.data[i][j][i3] << "," << simulation.data[i][j][i4] << ","
                << simulation.data[i][j][i5] << "," << simulation.data[i][j][i6] << endl;
        }
    }*/

    /*clog << simulation.data[0][0][simulation.descriptionIndex["node_angle_mean"]] << endl;
    simulation.calculateAllExtraData();
    clog << simulation.data[0][0][simulation.descriptionIndex["node_angle_mean"]] << endl;
    */

    simulation.deallocateData();
	return 0;
}

