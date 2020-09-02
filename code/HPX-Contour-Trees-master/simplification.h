#pragma once

#include <iostream>
#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/students_t.hpp>
#include <boost/math/distributions/chi_squared.hpp>
#include <boost/container/flat_set.hpp>

#include <random>
#include <set>

struct bqid {
    uint32_t value;
    uint32_t level;

    //most significant bits from 31 down to "level" are valid (inner nodes in the tree)
    //valid bits are compared, then a valid 1 is larger and a valid 0 is smaller than an invalid bit
    inline friend bool operator< (const bqid& lhs, const bqid& rhs){

        if (lhs.level == rhs.level){
            return lhs.value < rhs.value;
        }
        if (lhs.level > rhs.level){
            uint32_t left_actual_value = lhs.value >> lhs.level;
            uint32_t right_actual_value = rhs.value >> lhs.level;
            if (left_actual_value != right_actual_value ){
                return lhs.value < rhs.value;
            } else {
                return rhs.value & (1 << (lhs.level-1));
            }
        } else {
            uint32_t left_actual_value = lhs.value >> rhs.level;
            uint32_t right_actual_value = rhs.value >> rhs.level;
            if (left_actual_value != right_actual_value){
                return lhs.value < rhs.value;
            } else {
                return !(lhs.value & 1 << (rhs.level-1));
            }
        }
    }
    inline friend bool operator> (const bqid& lhs, const bqid& rhs){ return rhs < lhs; }
    inline friend bool operator<=(const bqid& lhs, const bqid& rhs){ return !(lhs > rhs); }
    inline friend bool operator>=(const bqid& lhs, const bqid& rhs){ return !(lhs < rhs); }
    inline friend bool operator==(const bqid& lhs, const bqid& rhs){ return !(lhs < rhs || rhs < lhs); }
    inline friend bool operator!=(const bqid& lhs, const bqid& rhs){ return !(lhs == rhs);}
};

struct bqnode {

    bqnode(){
        parent = nullptr;
        left = nullptr;
        right = nullptr;
        count = 0.0;
        L = 0.0;
    }

    bqnode* parent;
    bqnode* left;
    bqnode* right;
    double weight;
    double L;
    double count;
    bqid value;
};

struct bqnode_compare {
    bool operator() (const bqnode* lhs, const bqnode* rhs) const {
        return lhs->value < rhs->value;
    }
};

typedef boost::container::flat_set<bqnode*, bqnode_compare>::iterator iterator;

class bqsummary {

    public:
    boost::container::flat_set<bqnode*, bqnode_compare> bql;
    boost::container::flat_set<bqnode*, bqnode_compare> bqt;

    std::uint64_t n = 0;
    double alpha = 0.5;

    double L_bql = 0;
    double totalweight = 0;

    bqnode* dummy = new bqnode();

    bqnode* node(bqid val){
        dummy->value = val;
        return dummy;
    }

    bqsummary(double alpha){
        this->alpha = alpha;
    }

    void insert(uint32_t value){

        //inserted values are leaves, thus level = 0
        bqid val;
        val.value = value;
        val.level = 0;

        //increment overall sample size
        n++;
        totalweight++;

        //if bqt still empty or inserted value smaller than biggest bql entry
        iterator lb = bql.lower_bound(node(val));
        if (lb != bql.end() || bqt.size() == 0){
            //insert value to bql or increment counter of bql entry for value
            if ((lb != bql.end()) && (*lb)->value == val){
                (*lb)->count += 1.0;
                L_bql += 1.0;
            } else {
                bqnode* next = new bqnode();
                next->parent = nullptr;
                next->left = next;
                next->right = next;
                next->count = 1.0;
                next->value = val;
                bql.insert(next);
                L_bql += 1.0;
            }
        } else {
            //insert to bqt
            //from least common ancestor of value and biggest  bql entry
            bqid u = (*bql.nth(bql.size()-1))->value;
            bqid current = right(lca(u, val));
            //if no bqt entry from this ancestor exists insert it
            if (!bqt_contains(current)){
                bqnode* next = new bqnode();
                next->parent = nullptr;
                next->right = nullptr;
                next->left = nullptr;
                next->value = current;
                next->count = 1.0;
                next->L = L_bql;
                bqt.insert(next);
            } else {
                //follow path from this ancestor to value until a node does not exist
                while (true){
                    bqid next;
                    if (val < current){
                        next = left(current);
                    } else if (current < val){
                        next = right(current);
                    } else {
                        //if value itself exists increment count
                        (*bqt.find(node(current)))->count++;
                        break;
                    }
                    if (!bqt_contains(next)){
                        //if least existing ancestor of value is full insert closer ancestor, else increment count
                        bqnode* target = *bqt.find((node(current)));
                        if (target->count + 1.0 > alpha*target->L){
                            bqnode* insert_next = new bqnode();
                            insert_next->parent = target;
                            insert_next->value = next;
                            if (val < target->value){
                                target->left = insert_next;
                            } else {
                                target->right = insert_next;
                            }
                            insert_next->count = 1.0;
                            insert_next->L = target->L;
                            bqt.insert(insert_next);
                            break;
                        } else {
                            target->count += 1.0;
                            break;
                        }
                    } else {
                        current = next;
                    }
                }
            }
        }
    }

    double rank(uint32_t value){
        bqid val;
        val.value = value;
        val.level = 0;

        if (bql.count(node(val))){
            double result = 0;
            iterator current = bql.find(node(val));
            for (iterator it = bql.begin(); (it != current); it++)
                result += (*it)->count;
            return result;
        }

        bqid u = (*bql.nth(bql.size()-1))->value;
        bqid current = right(lca(u, val));
        double A = 0;
        while (true){
            bqid next;
            if (val < current){
                next = left(current);
            } else if (current < val){
                next = right(current);
            } else {
                //if value itself exists increment count
                return (*bqt.find(node(current)))->L;
            }
            if (!bqt_contains(next)){
               return (*bqt.find(node(current)))->L;
            } else {
                A += (*bqt.find(node(current)))->count;
                current = next;
            }
        }
    }

    bool rank_enough(uint32_t value, double p){
        double valrank = rank(value);
        return valrank >= n*p;
    }

    void update(){

        //identify z and u
        if (bql.size() > 1){
            double minBQL = 1.0/alpha;
            double remainBQL = 0.0;

            iterator cut;
            for (cut = bql.begin(); (cut != bql.end()); cut++){
                remainBQL += (*cut)->count;
                if (remainBQL >= minBQL)
                    break;
            }

            if (cut == bql.end())
                return;
            bqid u = (*cut)->value;
            cut++;
            if (cut == bql.end())
                return;
            bqid z = (*bql.nth(bql.size()-1))->value;

            //insert missing parents (all on path from z to lca(z,u)
            bqid rot = lca(z, u);
            bqid current;
            for (current = parent(z); (current != rot); current = parent(current)){
                if (bqt_contains(right(current))){
                    bqnode* dangling = *bqt.find(node(right(current)));
                    bqnode* next = new bqnode();
                    next->parent = nullptr;
                    next->left = nullptr;
                    next->right = dangling;
                    dangling->parent = next;
                    next->count = 0.0;
                    next->value = current;
                    next->L = 0.0;
                    bqt.insert(next);
                }
            }

            //insert elements from u through z from bql to bqt
            for (iterator it = cut; (it != bql.end()); ++it){
                L_bql -= (*it)->count;
                current = right(lca(u, (*it)->value));
                if (!(bqt_contains(current))){
                    bqnode* next = new bqnode();
                    next->parent = nullptr;
                    next->right = nullptr;
                    next->left = nullptr;
                    next->value = current;
                    next->count = 1.0;
                    next->L = L_bql;
                    bqt.insert(next);
                } else {
                    while (true){
                        bqid next;
                        if ((*it)->value < current){
                            next = left(current);
                        } else if (current < (*it)->value) {
                            next = right(current);
                        } else {
                            //leaf already exists
                            //should never happen here, we insert from bql thus we are strictly smaller than anything in bqt
                            (*bqt.find(node(current)))->count += 1.0;
                            break;
                        }
                        if (!bqt_contains(next)){
                            bqnode* target = *bqt.find(node(current));
                            if (target->count + 1.0 > alpha*target->L){
                                bqnode* insert_next = new bqnode();
                                insert_next->parent = target;
                                if ((*it)->value < target->value){
                                    insert_next->value = left(target->value);
                                    target->left = insert_next;
                                } else {
                                    insert_next->value = right(target->value);
                                    target->right = insert_next;
                                }
                                insert_next->count = 1.0;
                                insert_next->L = target->L;
                                bqt.insert(insert_next);
                                break;
                            } else {
                                target->count += 1.0;
                                break;
                            }
                        } else {
                            current = next;
                        }
                    }
                }
            }
            bql.erase(cut, bql.end());
        }

        calcweights();

        bqid u = (*bql.nth(bql.size()-1))->value;
        double L = L_bql;

        for (bqid current = parent(u); (current != parent(current)); current = parent(current)){
            if (bqt_contains(right(current))){
                compress(right(current),0, L);
                L += (*bqt.find(node(right(current))))->weight;
            }
        }
    }

    double totalcount(){

    }

    bool count(){
        double count = 0;

        for (iterator it = bql.begin(); (it != bql.end()); it++){
            bqnode* current = (*it);
            count += current->count;
        }

        assert(count == L_bql);

        for (iterator it = bqt.begin(); (it != bqt.end()); it++){
            bqnode* current = (*it);
            count += current->count;
        }

        assert(std::abs(count-n) < 0.01*n);

        if (bqt.size() > 0){
            iterator last = bqt.nth(bqt.size()-1);
            return((*last)->L + (*last)->weight == n);
        }
        return true;
    }

    void print(){

        std::cout << "N: " << n << std::endl;

        std::cout << "BQL Size: " << bql.size() << " BQL_L: " << L_bql << " (value, level, count, L)" << std::endl;

        for (iterator it = bql.begin(); (it != bql.end()); it++){
            bqnode* current = (*it);
            std::cout << "(" << current->value.value << ", " << current->value.level << ", " << current->count << ", " << current->L << ") ";
        }

        std::cout << std::endl << "BQT Size: " << bqt.size() << "(value, level, count, L)" << std::endl;

        for (iterator it = bqt.begin(); (it != bqt.end()); it++){
            bqnode* current = (*it);
            std::cout << "(" << current->value.value << ", " << current->value.level << ", " << current->weight << ", " << current->L << ") ";
        }
    }

    void calcweights(){
        bqid u = (*bql.nth(bql.size()-1))->value;
        bqid current;
        totalweight = L_bql;
        for (current = parent(u); (current != parent(current)); current = parent(current)){
            if (bqt_contains(right(current))){
                totalweight += weight(right(current));
            }
        }
    }

    double weight(bqid target){
        if (bqt_contains(target)){
            bqnode* current = *bqt.find(node(target));
            if (left(target) == right(target)){
                current->weight = current->count;
            } else {
                current->weight = current->count + weight(left(target)) + weight(right(target));
            }
            return current->weight;
        } else {
            return 0.0;
        }
    }

    void compress(bqid startnode, double dbt, double L){
        bqnode* current = *bqt.find(node(startnode));


        current->L = L;

        double c = current->count;

        if (left(startnode) == startnode){
            current->count -= dbt;
        } else {
            current->count = std::min(alpha * L, current->weight - dbt);
            if (current->count >= alpha*L) {
                dbt = dbt + current->count - c;
                double wl = 0.0;
                if (bqt_contains(left(startnode))){
                    wl = (*bqt.find(node(left(startnode))))->weight;
                    compress(left(startnode), std::min(dbt, wl), L);
                }
                if (bqt_contains(right(startnode)))
                    compress(right(startnode), std::max(dbt-wl, 0.0), L+wl+c+std::max(dbt-wl, 0.0));
            } else {
                erase(left(startnode));
                erase(right(startnode));
            }
        }
        if (current->count < 0.000001){
            erase(startnode);
        }
    }

    void erase(bqid target){
        if (!bqt_contains(target))
            return;
        if (left(target) != right(target)){
            erase(left(target));
            erase(right(target));
        }
        iterator it = bqt.find(node(target));
        bqnode* current = *it;
        bqt.erase(it);
        delete current;
    }

    bqid parent(bqid node){
        if (node.level == 32)
            return node;
        bqid result;
        result.value = node.value & (~(1 << node.level));
        result.level = node.level+1;
        return result;
    }

    bqid right(bqid node){
        if (node.level == 0)
            return node;

        bqid result;
        result.level = node.level - 1;
        result.value = node.value | (1 << result.level);
        return result;
    }

    bqid left(bqid node){
        if (node.level == 0)
            return node;

        bqid result;
        result.level = node.level - 1;
        result.value = node.value;
        return result;
    }

    bool bqt_contains(bqid target){
        iterator tmp = bqt.find(node(target));
        return (tmp != bqt.end());
    }

    bqid lca(bqid maxbql, bqid value){
        uint32_t common = ~(maxbql.value^value.value);
        uint32_t level;
        uint32_t resultval = 0;
        for (level = 31; (level >= 0); level--){
            if (!(common & (1 << level))){
                break;
            } else {
                resultval |= (maxbql.value & (1 << level));
            }
        }
        bqid result;
        result.value = resultval;
        result.level = level+1;
        return result;
    }
};

class Simplifier {

public:

    double significance = 0.01;
    double k = 3;
    double p = 0.95;

    double empiric_mean = 0;
    double empiric_variance = 0;
    std::uint64_t sample_size = 0;

    double T;
    double w = 0;
    double lower = std::numeric_limits<double>::max();
    double upper = std::numeric_limits<double>::max();

    int histogram_buckets;
    double histogram_min;
    double histogram_max;
    double bucket_width;

    std::vector<int> histogram;
    std::vector<int> ahistogram;
    std::vector<double> ahistogramValues;

    std::default_random_engine* generator;
    std::normal_distribution<float>* distribution;

    int counter = 0;

    bqsummary bq;

    Simplifier() : bq(0.002){

//        students_t t(n-1);
//        double T = quantile(complement(t, a/4));
//        double w = T * s/sqrt(double(n));

//        chi_squared chi(n-1);
//        double lower = sqrt((n-1) * s*s / quantile(complement(chi, a/4)));
//        double upper = sqrt((n-1)* s*s / quantile(chi, a/4));
    }

    void init_histogram(){
        histogram.resize(histogram_buckets, 0);
        bucket_width = (histogram_max - histogram_min)/(double)histogram_buckets;

        generator = new std::default_random_engine();
        distribution = new std::normal_distribution<float>(5.0,2.0);

        ahistogram.resize(histogram_buckets, 0);
        ahistogramValues.resize(histogram_buckets, 0.0);
        std::cout << "bucket width: " << bucket_width;
    }

    void print_histogram(){
        /*
        std::cout << "Count: " << sample_size;
        for (int i = 0; (i < histogram_buckets); i++){
            std::cout << "Histogram[" << i << "] = " << histogram[i] << std::endl;
        }
        */
    }

    void print_ahistogram(){
        std::cout << "Count: " << sample_size;
        for (int i = 0; (i < histogram_buckets); i++){
            std::cout << "Histogram[" << i << "] = " << ahistogram[i] << std::endl;
        }
    }

    bool ahistogram_large_enough(double arcLength){
        int myBucket = histogram_buckets - 1;
        for (int i = 0; (i < histogram_buckets-1); i++){
            if (ahistogramValues[i] <= arcLength){
                double lowdist = arcLength - ahistogramValues[i];
                if (ahistogramValues[i] == 0.0)
                    lowdist = 0;
                double highdist = arcLength - ahistogramValues[i+1];
                if (ahistogramValues[i+1] == 0.0)
                    highdist = 0;
                if (highdist > lowdist)
                    myBucket = i+1;
                else
                    myBucket = i;
                break;
            }
        }

        ahistogram[myBucket]++;
        ahistogramValues[myBucket] += (arcLength - ahistogramValues[myBucket])/(1.0*ahistogram[myBucket]);
        sample_size++;
        int quantileBucket = histogram_buckets;
        int sum = 0;
        while ((double)sum < (double)sample_size *(double)(1-p)){
            quantileBucket--;
            sum += ahistogram.at(quantileBucket);
        }
        return myBucket >= quantileBucket;
    }

    bool histogram_large_enough(double arcLength){
        double correctedArcLength = arcLength*100000;
        int myBucket = std::floor((correctedArcLength - histogram_min)/bucket_width);
        if (myBucket >= histogram_buckets)
            myBucket = histogram_buckets-1;
        histogram.at(myBucket)++;
        sample_size++;
        int quantileBucket = histogram_buckets;
        int sum = 0;
        while ((double)sum < (double)sample_size *(double)(1-p)){
            quantileBucket--;
            sum += histogram.at(quantileBucket);
        }
        return myBucket >= quantileBucket;

    }

    bool bqsummary_large_enough(float arcLength){

        arcLength = distribution->operator () ((*generator));

        //return true;
        std::uint32_t normalizedArcLength = *(std::uint32_t*) &arcLength;
        counter++;
        //std::cout << "Double: " << arcLength << " Int: " << normalizedArcLength << std::endl;

        bq.insert(normalizedArcLength);
        if (counter > 20){
            bq.update();
            counter = 0;
        }

        //if (bq.n < 10000){
        //    return true;
        //}
        return bq.rank_enough(normalizedArcLength, p);
    }

    bool gaussian_large_enough(double arcLength){

        arcLength = distribution->operator ()((*generator));

        //Tukey's Fences
        //return true;
        boost::math::normal_distribution<double> norm(empiric_mean-w, upper);
        double q1 = boost::math::quantile(norm, 0.25);
        double q3 = boost::math::quantile(norm, 0.75);
        double lower_fence = q1 - k*(q3-q1);
        double upper_fence = q3 + k*(q3-q1);

        if (arcLength < lower_fence)
            return false;

        if (arcLength > upper_fence)
            return true;

        //Numerically Stabilized Steiner Translation
        sample_size++;
        double diff = arcLength - empiric_mean;
        empiric_mean += (diff)/(1.0*sample_size);
        diff = arcLength - empiric_mean;
        empiric_variance += (diff)*(diff/(1.0*sample_size))*(sample_size-1.0);


        if (sample_size > 2){
            //Students t mean estimation
            boost::math::students_t t(sample_size-1);
            T = boost::math::quantile(boost::math::complement(t, significance/2));
            w = std::sqrt(T*T*empiric_variance/((sample_size-1)*sample_size));

            //Chi Squared standard-deviation estimation
            boost::math::chi_squared chi(sample_size-1);
            lower = std::sqrt(empiric_variance / boost::math::quantile(boost::math::complement(chi, significance/2)));
            upper = std::sqrt(empiric_variance / boost::math::quantile(chi, significance/2));

            //Estimated Normal Distribution Quantile Test
            boost::math::normal_distribution<double> norm(empiric_mean-w, upper);
            return arcLength > boost::math::quantile(norm, p);

        } else
            return false;
    }

};
