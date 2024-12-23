#include <unordered_map>
#include <fstream>
#include <iostream>
#include <sstream>
#include <chrono>
#include <vector>
#include <algorithm>
#include <queue>
#include <cmath>

using namespace std;

int NUM_OF_RESULTS = 10;
double BM25_K = 1.5;
double BM25_B = 0.75;

//struct definitions
struct LexiconEntry {
    uint64_t offset;     // Byte offset in the final index file
    uint32_t length;     // Number of bytes for this term's entry
    uint32_t docFreq;    // Number of documents containing the term
};

struct BlockMetaData {
    uint64_t offset;    // start position of the block in the index file
    uint32_t length;    // number of bytes of this block
    int lastDocID;  //lastdocid on this block
};

struct QueryResult {
    int docID;
    double score;
};

struct Posting {
    string term;
    int docID;
    double score;
};

struct Page {
    int termNum;
    uint64_t offset;
    uint64_t length;
};

//comparator for QueryResult min-heap based on score
struct CompareQueryScore {
    bool operator()(const QueryResult& a, const QueryResult& b) {
        return a.score > b.score; //lower score is on the top
    }
};

//comparator struct for sorting inverted index list by the docFreq in lexicon of this term
struct ListLengthComparator {
    const unordered_map<string, LexiconEntry>& lexicon;

    //constructor to initialize the lexicon reference
    ListLengthComparator(const unordered_map<string, LexiconEntry>& lex) : lexicon(lex) {}

    //comparison operator
    bool operator()(const pair<string, vector<uint8_t>>& a,
                    const pair<string, vector<uint8_t>>& b) const {
        //use the lexicon to get the length of list (in postings) for comparison
        return lexicon.at(a.first).docFreq < lexicon.at(b.first).docFreq;
    }
};

struct ComparePosting {
    bool operator()(const Posting &a, const Posting &b) {
        return a.docID > b.docID; // Min-heap based on score
    }
};


class TopResults {
private:
    priority_queue<QueryResult, vector<QueryResult>, CompareQueryScore> minHeap;
    int size;
public:
    TopResults(int k) {
        size = k;
    }

    //add to the min heap
    void add(const QueryResult& result) {
        //insert if less than 10 results
        if (minHeap.size() < size) {
            minHeap.push(result);
        } 
        //replace the lowest score element if there are already 10 elements
        else if (result.score > minHeap.top().score) {
            minHeap.pop(); 
            minHeap.push(result);
        }
    }

    //export min heap to vector
    vector<QueryResult> getResults() {
        vector<QueryResult> results;
        priority_queue<QueryResult, vector<QueryResult>, CompareQueryScore> intermediateHeap = minHeap;

        while (!intermediateHeap.empty())
        {
            results.push_back(intermediateHeap.top());
            intermediateHeap.pop();
        }

        reverse(results.begin(), results.end());
        return results;
    }

    double getLowestScore() {
        return minHeap.top().score;
    }

    bool isEmpty() {
        return minHeap.empty();
    }
};

unordered_map<string, LexiconEntry> loadLexicon(const string& filePath) {
    unordered_map<string, LexiconEntry> lexiconMap;
    ifstream lexiconFile(filePath);

    if (!lexiconFile.is_open()) {
        throw runtime_error("Failed to open lexicon file for reading: " + filePath);
    }

    string line;
    while (getline(lexiconFile, line)) {
        istringstream iss(line);
        string term;
        uint64_t offset;
        uint32_t length; 
        uint32_t docFreq;

        if (iss >> term >> offset >> length >> docFreq) {
            lexiconMap[term] = {offset, length, docFreq};
        } else {
            cerr << "Error parsing line: " << line << endl;
        }
    }

    lexiconFile.close();
    cout << "Lexicon file loaded" << endl;
    return lexiconMap;
}

unordered_map<int, Page> loadPageTable(const string& filePath) {
    unordered_map<int, Page> pageMap;
    ifstream pageTableFile(filePath);

    if (!pageTableFile.is_open()) {
        throw runtime_error("Failed to open page table file for reading: " + filePath);
    }

    string line;
    while (getline(pageTableFile, line)) {
        istringstream iss(line);
        int docID;
        int termNum; //length in words
        uint64_t offset;
        uint64_t length;

        if (iss >> docID >> termNum >> offset >> length) {
            Page currPage;
            currPage.termNum = termNum;
            currPage.offset = offset;
            currPage.length = length;
            pageMap[docID] = currPage;
        } else {
            cerr << "Error parsing line in page table file: " << line << endl;
        }
    }

    pageTableFile.close();
    cout << "Page table file loaded" << endl;
    return pageMap;
}

vector<BlockMetaData> loadBlockMetaData(const string& filePath) {
    vector<BlockMetaData> blockMetaDataVec;
    ifstream metaDataFile(filePath);

    if (!metaDataFile.is_open()) {
        throw runtime_error("Failed to open block metadata file for reading: " + filePath);
    }

    uint32_t offset = 0;
    string line;
    while (getline(metaDataFile, line)) {
        istringstream iss(line);
        uint32_t length;
        int lastDocID;

        if (iss >> length >> lastDocID) {
            blockMetaDataVec.push_back({offset, length, lastDocID});
            offset += length;
        } else {
            cerr << "Error parsing line in block meta data file: " << line << endl;
        }
    }

    metaDataFile.close();
    cout << "Block meta data file loaded" << endl;
    return blockMetaDataVec;
}

vector<uint8_t> openList(const string& term, const unordered_map<string, LexiconEntry> lexicon, ifstream& indexFile) {
    auto it = lexicon.find(term);
    //check if the term exists in the lexicon
    if (it == lexicon.end()) {
        cerr << "Term not found: " << term << endl;
        return {};
    }

    const LexiconEntry& entry = it->second;
    // seek to the start position for the term
    indexFile.seekg(entry.offset);
    // Read the specified number of bytes for the inverted index
    vector<uint8_t> buffer(entry.length);
    indexFile.read(reinterpret_cast<char*>(buffer.data()), entry.length);

    return buffer;

}

void sortListByLength(vector<pair<string, vector<uint8_t>>>& invertedLists, const unordered_map<string, LexiconEntry>& lexicon) {
    sort(invertedLists.begin(), invertedLists.end(), ListLengthComparator(lexicon));
}

//function to turn bytes to int
uint64_t byteToInt(const vector<uint8_t>& bytes) {
    uint64_t value = 0;
    int shift = 0;

    for (size_t i = 0; i < bytes.size(); ++i) {
        value |= (bytes[i] & 0x7F) << shift; // Mask the highest bit and shift
        shift += 7;

        // If the highest bit is not set, we are done
        if (!(bytes[i] & 0x80)) {
            break;
        }
    }

    return value;
}

//transform bytes to int vector
vector<int> bytesToIntVec(const vector<uint8_t>& bytes) {
    vector<int> numbers;
    vector<uint8_t> buffer;

    for (size_t i = 0; i < bytes.size(); ++i) {
        uint8_t byte = bytes[i];
        buffer.push_back(byte);

        // If the highest bit is not set, we have reached the end of the varbyte
        if (!(byte & 0x80)) {
            numbers.push_back(byteToInt(buffer)); // Convert to int and store
            buffer.clear(); // Clear buffer for the next varbyte
        }
    }
    
    // Handle case where buffer still has bytes (last varbyte)
    if (!buffer.empty()) {
        numbers.push_back(byteToInt(buffer));
    }

    return numbers;
}

//open lists for all the query terms
vector<pair<string, vector<uint8_t>>> readInvertedIndices(const vector<string>& terms, const unordered_map<string, LexiconEntry>& lexicon, const string& filePath) {
    vector<pair<string, vector<uint8_t>>> invertedLists;
    ifstream indexFile(filePath, ios::binary);

    if (!indexFile.is_open()) {
        throw runtime_error("Failed to open inverted index file for reading: " + filePath);
    }

    for (const string& term : terms) {
        vector<uint8_t> list = openList(term, lexicon, indexFile);
        if (!list.empty()) {
            pair<string,vector<uint8_t>> invertedList(term, list);
            invertedLists.push_back(invertedList);
        }
    }
    indexFile.close();
    return invertedLists;

}

//search the index of the block needed so that we can locate the block
int searchBlockIndex(const vector<BlockMetaData>& blockMetaDataVec, const int& listStartPos) {
    int left = 0;
    int right = blockMetaDataVec.size() -1;
    int index = -1;

    while (left <=right) {
        int mid = (left + right) /2;
        if (blockMetaDataVec[mid].offset < listStartPos) {
            left = mid + 1;
        }
        else if (blockMetaDataVec[mid].offset > listStartPos){
            right = mid - 1;
        }
        else {
            index = mid;
            break;
        }
    }
    return index;
}

template <typename T>
vector<T> sliceVector(const vector<T>& vec, int start, int length) {
    return vector<T>(vec.begin() + start, vec.begin()+start+length);
}

//binary search the next docID
int searchNextDocID(vector<int>docIDListBlock, int lookUpDocID) {
    int left = 0;
    int right = docIDListBlock.size() - 1;
    int foundIndex = -1;

    while (left <= right) {
        int mid = (left + right)/2;
        if (docIDListBlock[mid] < lookUpDocID) {
            left = mid + 1;
        }
        else {
            foundIndex = mid;
            right = mid - 1;
        }
    }
    return foundIndex;

}

void closeLists(vector<pair<string, vector<uint8_t>>>& invertedLists) {
    for (auto& invertedList : invertedLists) {
        (invertedList.second).clear();
    }
    invertedLists.clear();
}

pair<int,int> nextGEQ(const pair<string, vector<uint8_t>>& invertedList, 
            const int& lookUpDocID, const vector<BlockMetaData>& blockMetaDataVec, 
            const unordered_map<string, LexiconEntry>& lexiconMap) {
    //cout << "starting nextGEQ" << endl;
    //cout << "lookUpDocID = " << lookUpDocID << endl;
    string term = invertedList.first;
    
    vector<uint8_t> list = invertedList.second;
    int listStartPos = lexiconMap.at(term).offset;
    int listRestLength = lexiconMap.at(term).length;
    
    int foundIndex = -1;
    int foundDocID;
    int foundFreq;

    int blockIndex = searchBlockIndex(blockMetaDataVec, listStartPos);
    
    if (blockIndex == -1) {
        cerr << "error finding the next docID" << endl;
        return make_pair(-1,-1);
    }
    listStartPos = 0;

    bool startOfList = true;
    while (listRestLength >= 0) {
        if (blockMetaDataVec[blockIndex].lastDocID < lookUpDocID && listRestLength > blockMetaDataVec[blockIndex].length) {
            startOfList = false;
            
            listRestLength -= blockMetaDataVec[blockIndex].length;
            listStartPos += blockMetaDataVec[blockIndex].length;
            blockIndex++;
        }
        else {
            vector<uint8_t> listBlock = sliceVector(list, listStartPos, blockMetaDataVec[blockIndex].length);
            vector<int> decompressedListBlock = bytesToIntVec(listBlock);
            
            int listBlockSize = decompressedListBlock.size();
            
            vector<int> docIDListBlock = sliceVector(decompressedListBlock, 0, listBlockSize/2);
            
            int prevtDocID;
            if (startOfList) {
                prevtDocID = 0;
            } else {
                prevtDocID = blockMetaDataVec[blockIndex-1].lastDocID;
            }
            for (auto& docID : docIDListBlock)  {
                docID += prevtDocID;
                prevtDocID = docID;
            }
            
            if (docIDListBlock[docIDListBlock.size()-1] < lookUpDocID) {
                //cout << "reached the end" << endl;
                foundDocID = -1;
                foundFreq = -1;
                break;
            } else {
                foundIndex = searchNextDocID(docIDListBlock, lookUpDocID);
                foundDocID = docIDListBlock[foundIndex];
                foundFreq = decompressedListBlock[foundIndex + listBlockSize/2];
                break;
            }
            
            
        }
    }
    if (foundIndex == -1) {
        return make_pair(-1,-1);
    }
    else {
        return make_pair(foundDocID, foundFreq);
    }
}

double getScore(const string& term, const pair<int, int> docIDFreq, int& totalNumDoc, int& aveDocLength, const unordered_map<string, LexiconEntry>& lexiconMap, unordered_map<int, Page>& pageMap) {
    int numDocForTerm = lexiconMap.at(term).docFreq;
    
    int termFreq = docIDFreq.second;
    double score;
    score = log((totalNumDoc - numDocForTerm + 0.5) / (numDocForTerm + 0.5) + 1);
    
    score *= termFreq * (BM25_K+1)/(termFreq + BM25_K * (1-BM25_B + BM25_B * pageMap.at(docIDFreq.first).termNum / aveDocLength));
    return score;
}

vector<Posting> decompressWholeInvertedList(const pair<string, vector<uint8_t>>& invertedList, 
    const vector<BlockMetaData>& blockMetaDataVec, 
    const unordered_map<string, LexiconEntry>& lexiconMap,
    int& totalNumDoc, 
    int& aveDocLength, 
    unordered_map<int, Page>& pageMap) {

    vector<Posting> decompressedFullList;
    string term = invertedList.first;
    vector<uint8_t> list = invertedList.second;
    int listStartPos = lexiconMap.at(term).offset;
    int listendPosNext = listStartPos + lexiconMap.at(term).length;

    int blockIndex = searchBlockIndex(blockMetaDataVec, listStartPos);

    if (blockIndex == -1) {
        cerr << "error finding the next docID" << endl;
        return {{"-1",-1,-1}};
    }

    while (blockMetaDataVec[blockIndex].offset < listendPosNext) {
    
        vector<uint8_t> listBlock = sliceVector(list, blockMetaDataVec[blockIndex].offset-listStartPos, blockMetaDataVec[blockIndex].length);
        vector<int> decompressedListBlock = bytesToIntVec(listBlock);

        int prevDocID;
        int listBlockSize = decompressedListBlock.size();

        if (blockMetaDataVec[blockIndex].offset == listStartPos) {
            prevDocID = 0;
        }
        else {
            prevDocID = blockMetaDataVec[blockIndex - 1].lastDocID;
        }
        
        for (int i = 0; i < listBlockSize/2; i++) {
            int currDocID = decompressedListBlock[i] + prevDocID;
            int currFreq = decompressedListBlock[i + listBlockSize/2];
            
            double currScore = getScore(term, make_pair(currDocID, currFreq), totalNumDoc, aveDocLength, lexiconMap, pageMap);
            Posting currPosting;
            currPosting.term = term;
            currPosting.docID = currDocID;
            currPosting.score = currScore;
            
            decompressedFullList.push_back(currPosting);
            prevDocID = currPosting.docID;
        }
        blockIndex++;
        if (blockIndex >= blockMetaDataVec.size()) {
            break;
        }
    }
    return decompressedFullList;
}

Posting nextGEQ(const string& term, const vector<Posting>& decompressedFullList, 
            const int& lookUpDocID) {
    int left = 0;
    int right = decompressedFullList.size() - 1;
    int foundIndex = -1;

    while (left <= right) {
        int mid = (left + right)/2;
        if (decompressedFullList[mid].docID < lookUpDocID) {
            left = mid + 1;
        }
        else {
            foundIndex = mid;
            right = mid - 1;
        }
    }

    if (foundIndex == -1) {
        return {"-1", -1, -1};
    }
    return {term, decompressedFullList[foundIndex].docID, decompressedFullList[foundIndex].score};
}

bool isValidTerm(const string& queryTerm) {
    for (const char& c : queryTerm) {
        if (!isalnum(c) || static_cast<unsigned char>(c) >127 ) {
            return false;
        }
    }
    return true;
}

string lower(string& queryTerm) {
    string lowerCaseQuerTerm = "";
    for (char c : queryTerm) {
        lowerCaseQuerTerm += tolower(c);
    }
    return lowerCaseQuerTerm;
}

void andQuery(vector<pair<string, vector<uint8_t>>>& invertedLists, 
            const unordered_map<string, LexiconEntry>& lexiconMap, 
            const vector<BlockMetaData>& blockMetaDataVec, 
            unordered_map<int, Page>& pageMap, int& totalNumDoc, int& aveDocLength,
            TopResults& topResults) {
    //AND query
    cout << "AND QUERY" << endl;
    sortListByLength(invertedLists, lexiconMap);

    int docID = -1;
    Posting found;
    vector<pair<int,double>> foundList;
    if (invertedLists.empty()) {
        cout << "no matching found" << endl;
    }
    else {
        vector<Posting> decompressedFullList = decompressWholeInvertedList(invertedLists[0], blockMetaDataVec, lexiconMap, totalNumDoc, aveDocLength, pageMap);
    
        while (docID != 0) {
            string firstTerm = invertedLists[0].first;
            
            
            found = nextGEQ(firstTerm, decompressedFullList, docID);
            
            docID = found.docID;
            if (docID == -1) {
                break;
            }
            double score = found.score;
            pair<int, double> foundScore = make_pair(docID, score);
            
            
            foundList.push_back(foundScore);
            for (int i = 1; i < invertedLists.size(); i++) {
                
                pair<int, int> foundNextTerm = nextGEQ(invertedLists[i], docID, blockMetaDataVec, lexiconMap);
                
                //no intersection
                if (foundNextTerm.first > docID) {
                    docID = foundNextTerm.first;
                    foundList.clear();
                    docID--;
                    //cout << "no intersection" << endl;
                    break;
                }
                else {
                    //there is an intersection
                    if (foundNextTerm.first == -1) {
                        break;
                    }
                    double score = getScore(invertedLists[i].first, foundNextTerm, totalNumDoc, aveDocLength, lexiconMap, pageMap);
                    pair<int, double> foundNextTermScore = make_pair(docID, score);
                    foundList.push_back(foundNextTermScore);
                    
                    //cout << "found " << foundNextTerm.first << " " << foundNextTerm.second << endl;
                }
                
                
            }
            

            if (!foundList.empty()) {
            double score = 0;
                for (const auto& entry : foundList) {
                    score += entry.second;
                }
            QueryResult result;
            result.docID = foundList[0].first;
            result.score = score;
            topResults.add(result);
            }
            foundList.clear();

            docID++;
        }
    }
}

void orQuery(vector<pair<string, vector<uint8_t>>>& invertedLists, 
            const unordered_map<string, LexiconEntry>& lexiconMap, 
            const vector<BlockMetaData>& blockMetaDataVec, 
            unordered_map<int, Page>& pageMap, int& totalNumDoc, int& aveDocLength,
            TopResults& topResults) {
    //OR query
    cout << "OR QUERY" << endl;
    vector<vector<Posting>> postings;
    for (const auto& termList : invertedLists) {
        vector<Posting> currPostings = decompressWholeInvertedList(termList, blockMetaDataVec, lexiconMap, totalNumDoc, aveDocLength, pageMap);
        postings.push_back(currPostings);
    }

    priority_queue<Posting, vector<Posting>, ComparePosting> orMinHeap;
    vector<int> indices(postings.size(), 0);
    vector<double> maxScores(postings.size());

    //get maxScore for each list
    //initialize minHeap for or query
    for (int i = 0; i < postings.size(); i++) {
        if (!postings[i].empty()) {
            maxScores[i] = postings[i][0].score;
            for (const auto &posting : postings[i]) {
                maxScores[i] = max(maxScores[i], posting.score); 
            }
            orMinHeap.push(postings[i][0]); //push the first entry of each postings list
        }
    }

    /*for (int i = 0; i < maxScores.size(); i++) {
        cout << "maxScore is " << maxScores[i] << endl;
    }*/

    while (!orMinHeap.empty()) {
        Posting currentPosting = orMinHeap.top();
        orMinHeap.pop();
        
        double totalScore = currentPosting.score;

        //add same docID's score
        while (!orMinHeap.empty() && orMinHeap.top().docID == currentPosting.docID) {
            
            totalScore += orMinHeap.top().score;
            orMinHeap.pop();    
        }

        QueryResult result;
        result.docID = currentPosting.docID;
        result.score = totalScore;
        topResults.add(result);

        //maxScore early termination, not really working right now
        //we don't really use maxScore now
        double lowestScoreInResults = topResults.getLowestScore(); //the smallest score in the results
        for (int i = 0; i < postings.size(); ++i) {
            if (indices[i] < postings[i].size() && postings[i][indices[i]].docID == currentPosting.docID) {
                indices[i]++; // Move index for this term
                if (indices[i] < postings[i].size()) {
                    orMinHeap.push(postings[i][indices[i]]);   
                }
            }
        }
    }
}


int main() {

    
    string lexiconFilePath = "output/lexicon.txt";
    string pageTableFilePath = "output/pagetable.tsv";
    string blockMetaDataFilePath = "output/blockMetaData.txt";
    string indexFilePath = "output/index.bin";
    string dataFilePath = "collection/collection.tsv";
    
    try {

        ifstream inFile(dataFilePath);

        if (!inFile.is_open()) {
            cout << "fail to read the dataset" << endl;
            return 0;
        }
        
        unordered_map<string, LexiconEntry> lexiconMap = loadLexicon(lexiconFilePath);

        
        unordered_map<int, Page> pageMap = loadPageTable(pageTableFilePath);
        int totalNumDoc = pageMap.size();
        int totalDocLength = 0;
        for (const auto& page : pageMap) {
            totalDocLength += page.second.termNum;
        }
        int aveDocLength = totalDocLength / totalNumDoc;

        
        vector<BlockMetaData> blockMetaDataVec = loadBlockMetaData(blockMetaDataFilePath);
    

        cout << "Search Engine is Ready" << endl;
        cout << endl;

        string userInput;
        bool inputValid = true;

        while(true) {
            TopResults topResults(NUM_OF_RESULTS);
            cout << "What do you want to search about?" << endl;
            cout << "Please separate your terms with space. The term needs to be consisted of digit or letter" << endl;
            cout << "Enter '-1' to exit" << endl;
            getline(cin, userInput);

            if (userInput == "-1") {
                cout << "Query done!" << endl;
                break;
            }

            istringstream inputStream(userInput);
            string queryTerm;
            vector<string> queryTerms;

            while (inputStream >> queryTerm) {
                if (isValidTerm(queryTerm)) {
                    queryTerms.push_back(lower(queryTerm));
                } else {
                    cout << "invalid input" << endl;
                    inputValid = false;
                    continue;
                }
            }

            if (!inputValid) {
                continue;
            }
            

            //read in inverted index lists
            vector<pair<string, vector<uint8_t>>> invertedLists = readInvertedIndices(queryTerms, lexiconMap, indexFilePath);
            cout << endl;

            while (userInput != "1" && userInput != "2") {
                cout << "AND (type 1) or OR (type 2) query?" << endl;
                getline(cin, userInput);
            }

            cout << endl;

            // Start the timer
            auto start = chrono::high_resolution_clock::now();
            if (userInput == "1") {
                andQuery(invertedLists, lexiconMap, blockMetaDataVec, pageMap, totalNumDoc, aveDocLength, topResults);
            } else {
                orQuery(invertedLists, lexiconMap, blockMetaDataVec, pageMap, totalNumDoc, aveDocLength, topResults);
            }

            vector<QueryResult> results = topResults.getResults();
            cout << "RESULTS ARE:" << endl;
            if (topResults.isEmpty()) {
                cout << "No Results" << endl;
            }
            else {
                for (const auto& result : results) {
                    cout << "DocID = " << result.docID <<'\t' << "Score = " << result.score << endl;
                    //read the corresponding line from dataset
                    uint32_t lineLength;
                    inFile.seekg(pageMap[result.docID].offset);
                    lineLength = pageMap[result.docID].length;
                    string line(lineLength, '\0');
                    inFile.read(&line[0], lineLength);
                    //extract the text
                    size_t tabPos = line.find('\t');
                    string text = "";
                    if (tabPos != string::npos) {  
                        text =  line.substr(tabPos + 1);
                    }
                    if (text.length() != 0) {
                        cout << text << endl;
                        cout << endl;
                    }
                }
            }
            cout << endl;

            // Stop the timer
            auto end = chrono::high_resolution_clock::now();
            chrono::duration<double> duration = end - start;
            cout << "Query time: " << duration.count() << " seconds." << endl;

            closeLists(invertedLists);

        }
    } catch (const exception& e) {
        cerr << e.what() << endl;
    }
    
    return 0;

}


