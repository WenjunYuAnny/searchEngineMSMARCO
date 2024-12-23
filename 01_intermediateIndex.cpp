// indexer.cpp
#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <vector>
#include <utility>
#include <cctype>
#include <cstdlib>
#include <unordered_set>
#include <unordered_map>
#include <filesystem> // Include filesystem library
#include <chrono>

namespace fs = std::filesystem;
using namespace std;

const int MAX_BLOCK_SIZE = 100 * 1024 * 1024; // 100 MB per block

// Posting structure: docID and term frequency
struct Posting {
    int docID;
    int termFreq;
};

//Page structure
struct Page {
    int docID;
    int termNum;
    uint64_t offset; //start position of the line in the file
    uint64_t length; // length of the line
};

// Function to tokenize text ASCII check
vector<std::string> tokenize(const string& text) {
    std::vector<std::string> tokens;
    std::string token;
    auto isASCII = [](const string& token) -> bool {
        for (char c : token) {
            if (static_cast<unsigned char>(c) > 127) return false;
        }
        return true;
    };

    for (char c : text) {
        if (std::isalnum(static_cast<unsigned char>(c))) {
            token += std::tolower(static_cast<unsigned char>(c));
        } else if (!token.empty()) {
            if (isASCII(token)) {
                tokens.push_back(token);
            }
            token.clear();
        }
    }
    if (!token.empty() && isASCII(token)) {
        tokens.push_back(token);
    }
    return tokens;
}

// Function to write intermediate posting file in text format
void writeTextPostingFile(const string& filename,
                          const map<string, vector<Posting>>& invertedIndex) {
    ofstream outfile(filename);
    if (!outfile.is_open()) {
        throw runtime_error("Failed to open text intermediate file for writing: " + filename);
    }

    for (const auto& [term, postings] : invertedIndex) {
        // Write the term
        outfile << term;

        // Write each posting as docID:termFreq separated by space
        for (const auto& posting : postings) {
            outfile << " " << posting.docID << ":" << posting.termFreq;
        }

        // End the line for the current term
        outfile << "\n";
    }
    outfile.close();
}

void writePageTableFile(const string& pageTableFileName, const vector<Page>& pageTableVec) {
    ofstream outfile(pageTableFileName);
    if (!outfile.is_open()) {
        throw runtime_error("Failed to open page table file for writing: " + pageTableFileName);
    }

    for (const auto& page : pageTableVec) {
        outfile << page.docID << '\t' << page.termNum << '\t' << page.offset << '\t' << page.length << endl;
    }

    outfile.close();
    cout << "Page Table File Completed" << endl;
}

// Function to parse the collection and create intermediate posting files
void parseCollectionWritePageTable(const string& inputFilePath, 
                    map<std::string, vector<Posting>>& invertedIndex,
                    int& currentBlockSize,
                    const int maxBlockSize,
                    int& blockCount,
                    const string& outputDir,
                    const string& pageTableFileName,
                    vector<Page> pageTableVec) {
    ifstream infile(inputFilePath);
    if (!infile.is_open()) {
        throw runtime_error("Failed to open collection file: " + inputFilePath);
    }

    string line;
    static int processedDocs = 0; // Document counter
    uint64_t startPos;
    while (getline(infile, line)) {
        // Split the line into docID and passage
        size_t tabPos = line.find('\t');
        if (tabPos == string::npos) {
            startPos += line.length() + 1;
            continue; // Skip malformed lines
        }
        int docID = stoi(line.substr(0, tabPos));
        string passage = line.substr(tabPos + 1);

        // Tokenize the passage
        vector<string> tokens = tokenize(passage);

        //wirte to page table file
        Page currPage;
        currPage.docID = docID;
        currPage.length = line.length();
        currPage.offset = startPos;
        currPage.termNum = tokens.size();
        pageTableVec.push_back(currPage);
        startPos += currPage.length + 1;
        
        // Count term frequencies in the current document
        unordered_map<string, int> termFreqMap;
        for (const auto& token : tokens) {
            termFreqMap[token]++;
        }

        // Update the inverted index
        for (const auto& [term, freq] : termFreqMap) {
            invertedIndex[term].emplace_back(Posting{docID, freq});
            // Estimate block size: term length + (docID + freq) bytes
            currentBlockSize += term.size() + sizeof(int) * 2;
        }

        // Increment document counter and log progress
        processedDocs++;
        if (processedDocs % 100000 == 0) {
            cout << "Processed " << processedDocs << " documents..." << endl;
        }

        // Check if the current block size exceeds the maximum allowed
        if (currentBlockSize >= maxBlockSize) {
            // Write the current inverted index to an intermediate file
            string filename;
            if (blockCount < 10) {
                filename = outputDir + "/intermediate_0" + to_string(blockCount++) + ".txt";
            } else {
                filename = outputDir + "/intermediate_" + to_string(blockCount++) + ".txt";
            }
            writeTextPostingFile(filename, invertedIndex);
            cout << "Written intermediate file: " << filename << endl;

            // Clear the in-memory inverted index and reset block size
            invertedIndex.clear();
            currentBlockSize = 0;
        }
    }

    infile.close();

    // Write any remaining postings to an intermediate file
    if (!invertedIndex.empty()) {
        string filename = outputDir + "/intermediate_" + to_string(blockCount++) + ".txt";
        writeTextPostingFile(filename, invertedIndex);
        cout << "Written intermediate file: " << filename << endl;
    }

    if (!pageTableVec.empty()) {
        writePageTableFile(pageTableFileName, pageTableVec);
    }       
} 

int main() {

    // Start the timer
    auto start = chrono::high_resolution_clock::now();

    string inputFilePath = "collection/collection.tsv";
    string outputDir = "temp";
    string pageTableFileName = "output/pagetable.tsv";

    // Check if output directory exists; if not, create it
    try {
        if (!fs::exists(outputDir)) {
            fs::create_directories(outputDir);
            std::cout << "Created output directory: " << outputDir << std::endl;
        } else {
            // Optional: Check if it's a directory
            if (!fs::is_directory(outputDir)) {
                std::cerr << "Output path exists and is not a directory: " << outputDir << std::endl;
                return EXIT_FAILURE;
            }
        }
    } catch (const fs::filesystem_error& e) {
        std::cerr << "Filesystem error: " << e.what() << std::endl;
        return EXIT_FAILURE;
    }

    // Initialize variables
    map<string, vector<Posting>> invertedIndex;
    vector<Page> pageTableVec;
    int currentBlockSize = 0;
    int blockCount = 0;

    try {
        parseCollectionWritePageTable(inputFilePath, invertedIndex, currentBlockSize, MAX_BLOCK_SIZE, blockCount, outputDir, pageTableFileName, pageTableVec);
        std::cout << "Indexing completed successfully. " << blockCount << " intermediate files created." << std::endl;
    } catch (const std::exception& ex) {
        std::cerr << "Error during indexing: " << ex.what() << std::endl;
        return EXIT_FAILURE;
    }

    // Stop the timer
    auto end = chrono::high_resolution_clock::now();
    chrono::duration<double> duration = end - start;
    cout << "Execution time: " << duration.count() << " seconds." << std::endl;

    return EXIT_SUCCESS;
}