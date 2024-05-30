//
// Created by polya on 5/26/24.
//

#ifndef T_PHASE_NO_GRAD_PARSEPKL_HPP
#define T_PHASE_NO_GRAD_PARSEPKL_HPP

#include "version1LJPotPBC2Atom.hpp"

//this subroutine parses xml or bin files, we parse bin files for speed
class reader {

public:
    reader(const int &rowNum,const int &TInd, const unsigned long long& cellNum) {
        this->cellNum=cellNum;
        this->rowNum = rowNum;
        this->TRoot = "./version1Data/1d/funcLJPotPBC/row" + std::to_string(rowNum) + "/";
        std::vector<std::string >sortedTDirs= scanTDirs(TRoot);
        this->TDir=TRoot+"/"+sortedTDirs[TInd]+"/";
        std::cout<<"selected TDir="<<TDir<<std::endl;
        try {
            this->UInOneFile = std::shared_ptr<double[]>(new double[version1dLJPot2Atom::loopMax],
                                                         std::default_delete<double[]>());
            this->x_inOneFile = std::shared_ptr<double[]>(new double[version1dLJPot2Atom::loopMax * cellNum],
                                                          std::default_delete<double[]>());
        }
        catch (const std::bad_alloc& e) {
            std::cerr << "Memory allocation error: " << e.what() << std::endl;
        } catch (const std::exception& e) {
            std::cerr << "Exception: " << e.what() << std::endl;
        }
    }

public:
    ///
    /// @param TPath path containing all T folders
    /// @return T folder names sorted by T
     std::vector<std::string> scanTDirs(const std::string &TPath){
        std::vector<std::string> TDirs;
        std::string searchPath = TPath;
        std::string TPattern = "T([+-]?(\\d+(\\.\\d*)?|\\.\\d+)([eE][+-]?\\d+)?)";
        std::vector<double> TValsAll;

        if (fs::exists(searchPath) && fs::is_directory(searchPath)) {
            for (const auto &entry : fs::directory_iterator(searchPath)) {
                if (entry.path().filename().string()[0] == 'T' and entry.path().filename().string()[1] != 'x') {
                    TDirs.push_back(entry.path().filename().string());
                    std::smatch matchT;
                    if (std::regex_search(entry.path().filename().string(), matchT, std::regex(TPattern))) {
                        double TVal = std::stod(matchT.str(1));
                        TValsAll.push_back(TVal);
                    }
                }
            }
        }

        std::vector<size_t> inds = argsort(TValsAll);
        std::vector<std::string> sortedFiles;
        for (const auto &i : inds) {
            sortedFiles.push_back(TDirs[i]);
        }

        return sortedFiles;


    }//end function scanTDirs



    template<class T>
    std::vector<size_t> argsort(const std::vector<T> &v) {
        std::vector<size_t> idx(v.size());
        std::iota(idx.begin(), idx.end(), 0);
        std::stable_sort(idx.begin(), idx.end(), [&v](size_t i1, size_t i2) { return v[i1] <= v[i2]; });
        return idx;
    }

    template<class T>
    static void printVec(const std::vector<T> &vec) {
        for (int i = 0; i < vec.size() - 1; i++) {
            std::cout << vec[i] << ",";
        }
        std::cout << vec[vec.size() - 1] << std::endl;
    }

    ///UAll, xA_All, xB_All folder's files
    void searchFiles();


    ///
    /// @param path the path containing xml files
    /// @return sorted bin files by starting loop or end loop
    std::vector<std::string> sortOneDir(const std::vector<std::string> &allFiles);

    ///sort files by starting loop
    void sortFiles();

    void parseSummary();

    ///
    /// @param filename file name
    /// @param values values in file
    /// @param number_of_values number of values
    bool loadMsgFile(const std::string& filename, std::shared_ptr<double[]>& values, size_t& number_of_values);


    unsigned long long loadU();

    unsigned long long load_xA();

    unsigned long long load_xB();

    ///data to json, json as input to plot
    void data2json();

    ///compute the column means of arma_xA, arma_xB
    void colmeans();

    static void printMat(const arma::dmat &mat, std::ostream & os)  {
        int n = mat.n_rows;
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n - 1; j++) {
                os << mat(i, j) << ",";
            }
            os<< mat(i, n - 1) << std::endl;
        }
    }

    ///compute correlation functions GAA
    void computeGAA();

    ///compute correlation functions GAB
    void computeGAB();

    ///compute correlation functions GBB
    void computeGBB();







public:
    int TInd;
    std::string TDir;
    std::string TRoot;
    int rowNum;
    unsigned long long cellNum;
    std::string UPath;
    std::string xAPath;
    std::string xBPath;

    unsigned long long nEqCounterStart=0;
    unsigned long long lagEst=0;

    std::vector<std::string> UFilesAll;
    std::vector<std::string> xAFilesAll;
    std::vector<std::string> xBFilesAll;
    std::vector<std::string> sorted_UFilesAll;
    std::vector<std::string> sorted_xAFilesAll;
    std::vector<std::string> sorted_xBFilesAll;

    std::shared_ptr<double[]> UInOneFile;
    std::shared_ptr<double[]>x_inOneFile;

    std::shared_ptr<double[]> USelected;
    std::shared_ptr<double[]>xA_selected;
    std::shared_ptr<double[]> xB_selected;
    arma::dmat arma_xA;
    arma::dmat arma_xB;
    arma::drowvec E_xARow;
    arma::dcolvec E_xACol;

    arma::drowvec E_xBRow;
    arma::dcolvec E_xBCol;

    arma::dmat E_xA2;
    arma::dmat E_xB2;
};

#endif //T_PHASE_NO_GRAD_PARSEPKL_HPP
