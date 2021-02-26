#pragma once

#include <vector>
#include <map>
#include <cmath>

//* Network class with merging cluster by Newman-Ziff algorithm
struct NZ_Network{
private:
    //* Size of nodes and links
    int m_size{0};
    int m_linkSize{0};

    //* Maximum Cluster and Second Maximum Cluster
    int m_maximumClusterSize{1};
    int m_secondMaximumClusterSize{1};
    int m_deltaMaximumClusterSize{0};

    //* m_parent[node] : parent of each 'node'
    std::vector<int> m_parent;

    //* m_sortedCluster[size] : number of cluster of 'size'
    std::map<int, int> m_sortedCluster;

    //* m_birth[root] : birth time of each 'root'
    //* changedAge[root] : {age, size} of cluster with 'root'
    std::vector<int> m_birth;
    std::vector<std::pair<int, int>> m_changedAge;

public:
    //* Constructor
    NZ_Network() {}

    NZ_Network(const int &t_size)
    : m_size(t_size)
    {
        //* Make every node to root node with size 1
        m_parent.resize(t_size,-1);

        //* Initialize Sorted Cluster
        m_sortedCluster[1] = t_size;

        //* Initialize birth time with 0 and changedAge
        m_birth.resize(t_size);
        m_changedAge.resize(2,std::pair<int, int> {0,0});
    }

    //* Simple get functions
    int getMaximumClusterSize() const {return m_maximumClusterSize;}
    int getClusterSize(const int& t_root) const {return -m_parent[t_root];}
    int getDeltaMaximumClusterSize() const {return m_deltaMaximumClusterSize;}
    std::vector<std::pair<int,int>> getChangedAge() const {return m_changedAge;}


    //* get the root of input node
    int getRoot(const int &t_node){
        //* t_node is root
        if (m_parent[t_node] < 0){
            return t_node;
        }
        //* recursively find node
        return m_parent[t_node] = getRoot(m_parent[t_node]);
    }

    //* Merge two clusters
    void merge(const int& t_root1, const int& t_root2){
        //! update link size
        m_linkSize++;

        //! Save age and reset birth time
        m_changedAge[0] = std::pair<int,int>{m_linkSize-m_birth[t_root1], -m_parent[t_root1]};
        m_changedAge[1] = std::pair<int,int>{m_linkSize-m_birth[t_root2], -m_parent[t_root2]};
        m_birth[t_root1] = m_linkSize;
        m_birth[t_root2] = m_linkSize;

        //! Get the size of each clusters
        const int size1 = -m_parent[t_root1];
        const int size2 = -m_parent[t_root2];
        const int newSize = size1+size2;

        //! Update Parent
        m_parent[t_root1] -= size2;
        m_parent[t_root2] = t_root1;

        //! Update Sorted Cluster
        --m_sortedCluster[size1];
        if (m_sortedCluster[size1]==0){
            m_sortedCluster.erase(size1);
        }
        --m_sortedCluster[size2];
        if (m_sortedCluster[size2]==0){
            m_sortedCluster.erase(size2);
        }
        ++m_sortedCluster[newSize];

        //! Find maximum cluster
        if (m_maximumClusterSize < newSize){
            m_deltaMaximumClusterSize = newSize-m_maximumClusterSize;
            m_maximumClusterSize = newSize;
        }
        else{
            m_deltaMaximumClusterSize = 0;
        }
    }

    //* Calculate second Maximum Cluster Size
    void processSecondMaximumClusterSize(){
        for (auto it=m_sortedCluster.rbegin(); it != m_sortedCluster.rend(); ++it){
            if (it->first != m_maximumClusterSize || it->second>1){
                m_secondMaximumClusterSize = it->first;
                break;
            }
        }
    }

    //* Get second Maximum Cluster Size
    int getSecondMaximumClusterSize(){
        processSecondMaximumClusterSize();
        return m_secondMaximumClusterSize;
    }

    //* Get sorted Cluster
    std::map<int,int> getSortedCluster(const int &excludeNum=1){
        std::map<int,int> result=m_sortedCluster;

        //! exclude maximum cluster
        --result[m_maximumClusterSize];

        //! exclude second giant
        if (excludeNum-1){
            processSecondMaximumClusterSize();
            --result[m_secondMaximumClusterSize];
        }
        return result;
    }

    //* Get Mean cluster size
    double getMeanClusterSize() const{
        const int firstMoment = m_size-m_maximumClusterSize;
        double secondMoment = 0;
        for (auto it=m_sortedCluster.begin(); it!=m_sortedCluster.end(); ++it){
            secondMoment += pow(it->first,2)*it->second;
        }
        //! exclude infinite size cluster
        secondMoment -= pow(m_maximumClusterSize,2);

        return secondMoment/firstMoment;
    }
};