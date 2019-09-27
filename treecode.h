//
// Created by lurker on 9/25/19.
//

#ifndef TREECODE_H
#define TREECODE_H

#if !defined __extern_always_inline && defined __clang__
# if defined __GNUC_STDC_INLINE__ || defined __GNUC_GNU_INLINE__
#  define __extern_inline extern __inline __attribute__ ((__gnu_inline__))
#  define __extern_always_inline \
  extern __always_inline __attribute__ ((__gnu_inline__))
# else
#  define __extern_inline extern __inline
#  define __extern_always_inline extern __always_inline
# endif
#endif

#include "utils.h"
#include "blas_wrapper.h"
#include "linalg.h"

class point {
public:
    scalar_t x;
    scalar_t y;

    point() : x(0.), y(0.) {}
    point(scalar_t _x, scalar_t _y) : x(_x), y(_y) {}

    virtual ~point() = default;

    bool operator>=(const point &a) {
        return (x >= a.x - EPS) && (y >= a.y - EPS);
    }

    bool operator<=(const point &a) {
        return (x <= a.x + EPS) && (y <= a.y + EPS);
    }

    bool operator==(const point &a) {
        return fabs(x - a.x) < EPS && fabs(y - a.y) < EPS;
    }

};

class baseNode {
public:
    index_t parent;
    index_t child[4]{};

    index_t nLevel;
    index_t nodeIndex;

    point center;
    point radius;

    index_t nSource;
    vector<index_t> sourceIndex;

    bool isLeaf;
    bool isEmpty;

    bool chargeComputed;

    baseNode(index_t level, index_t index) {
        parent = -1;
        for (int & i : child) {
            i = -1;
        }
        nLevel = level;
        nodeIndex = index;
        isLeaf = false;
        isEmpty = false;
        chargeComputed = false;
        nSource = 0;
    }

    virtual ~baseNode() = default;
};

class node : public baseNode {
public:
    node(index_t level, index_t index) : baseNode(level, index) {}
    ~node() override = default;

    vector<point> scaledCnode;
    Vector potential;
    Vector nodePotential;
    Vector charge;
    Vector nodeCharge;

    Matrix R;
};

class tree {
public:
    vector<node> dict;
    index_t maxId;
    index_t root;
    index_t nSource;

    index_t rank;
    index_t maxLevel;

    vector<point> sourceTree;

    point center;
    point radius;

    tree() {
        maxId = -1;
        root = -1;
        nSource = 0;
        rank = 0;
        maxLevel = 0;
    }
    ~tree() = default;

    void populate(vector<point> &_source, index_t _nSource, index_t _rank, index_t _maxLevel) {
        sourceTree = _source;
        nSource = _nSource;
        maxLevel = 0;
        rank = _rank;

        getCenterRadius(_source);

        root = 0;
        dict.clear();
        dict.emplace_back(0, 0);
        maxId = root;

        dict[root].nSource = nSource;
        dict[root].center = center;
        dict[root].radius = radius;
        dict[root].sourceIndex.resize((unsigned long) nSource);

        for (index_t i = 0; i < nSource; ++i) {
            dict[root].sourceIndex[i] = i;
        }
        RUN("initialization", assignChildren(root, _maxLevel));
    }

protected:
    void getCenterRadius(vector<point> &_source) {
        assert(_source.size() > 0);
        scalar_t x_max = _source[0].x;
        scalar_t x_min = _source[0].x;
        scalar_t y_max = _source[0].y;
        scalar_t y_min = _source[0].y;

        for (size_t i = 0; i < _source.size(); ++i) {
            x_max = std::max(x_max, _source[i].x);
            y_max = std::max(y_max, _source[i].y);
            x_min = std::min(x_min, _source[i].x);
            y_min = std::min(y_min, _source[i].y);
        }
        this->center.x = (x_max + x_min) / 2.0;
        this->center.y = (y_max + y_min) / 2.0;
        this->radius.x = (x_max - x_min) / 2.0;
        this->radius.y = (y_max - y_min) / 2.0;
    }

    void assignChildren(index_t _id, index_t _maxLevel) {
        /*
         * when assigning children nodes, the points are not assigned due to storage.
         *
         * Now the limitation of nodes is around 2^24.
         */
        assert(root != -1); // check tree is non-empty

        // check source
        if (dict[_id].nSource == 0) {
            dict[_id].isLeaf = true;
            dict[_id].isEmpty = true;
        } else {
            // divide
            if ((dict[_id].nSource <= rank) || (dict[_id].nLevel == _maxLevel)) {
                dict[_id].isLeaf = true;
                if (maxLevel < dict[_id].nLevel) {
                    maxLevel = dict[_id].nLevel;
                }
            } else {
                // not a leaf
                for (index_t i = 0; i < 4; ++i) {
                    maxId += 1;
                    dict[_id].child[i] = maxId;
                    dict.emplace_back(dict[_id].nLevel + 1, i);
                    dict[maxId].parent = _id;
                    dict[maxId].center.x = dict[_id].center.x + ((i & 1) - 0.5) * dict[_id].radius.x;
                    dict[maxId].center.y = dict[_id].center.y + (((i >> 1) & 1) - 0.5) * dict[_id].radius.y;
                    dict[maxId].radius.x = dict[_id].radius.x * 0.5;
                    dict[maxId].radius.y = dict[_id].radius.y * 0.5;
                    dict[maxId].nSource = 0;

                }

                /*
                 * can be accelerated by **reduce**
                 */
                for (index_t i = 0; i < dict[_id].nSource; ++i) {
                    index_t index = dict[_id].sourceIndex[i];
                    index_t y_bit = sourceTree[index].y < dict[_id].center.y ? 0 : 1;
                    index_t x_bit = sourceTree[index].x < dict[_id].center.x ? 0 : 1;
                    index_t childIndex = 2 * y_bit + x_bit;

                    index_t childId = dict[_id].child[childIndex];
                    dict[childId].sourceIndex.push_back(index);
                    dict[childId].nSource += 1;
                }

                for (index_t i = 0; i < 4; ++i) {
                    assignChildren(dict[_id].child[i], _maxLevel);
                }
            }
        }
    }
};

class TreeCode{

public:
    tree t;
    Vector chargeTree;

    std::function<scalar_t (point&, point&) > eval;
    index_t rank;

    bool selfExclusive;

    Matrix R[4];

    index_t nChebyshev;
    Vector chebyNode;
    Matrix tNode;

    TreeCode() {
        nChebyshev = 0;
        rank = 0;
        selfExclusive = false;
    }

    ~TreeCode() = default;

    void initalize(index_t _nChebyshev, vector<point> &_source, Vector &_charge,
            index_t _nSource, index_t  _rank, index_t _maxLevel) {

        t.populate(_source, _nSource, _rank, _maxLevel);

        nChebyshev = _nChebyshev;
        chargeTree = _charge;

        rank = nChebyshev * nChebyshev;

        chebyNode = Vector(nChebyshev);
        getStandardChebyNodes(_nChebyshev, chebyNode);

        tNode = Matrix(nChebyshev, nChebyshev);
        getStandardChebyPoly(nChebyshev, nChebyshev, chebyNode, tNode);



        getTransfer(nChebyshev, chebyNode, tNode, R);



    }

    void getStandardChebyNodes(index_t _nChebyshev, Vector &_chebyNode) {
        assert(_chebyNode.row() == nChebyshev);
        for (index_t i = 0; i < _nChebyshev; ++i) {
            _chebyNode(i) = -cos((i + 0.5) * M_PI / _nChebyshev);
        }
    }

    void getStandardChebyPoly(index_t _nChebyPoly, index_t _N, Vector &_x, Matrix &_T) {
        assert(_T.row() == _N);
        assert(_T.col() == _nChebyPoly);

        setValue(_T, 0);

        Vector ones(_N);
        setValue(ones, 1.0);
        _T.setColumn(0, ones);

        if (_nChebyPoly > 1) {
            _T.setColumn(1, _x);
            for (index_t i = 2; i < _nChebyPoly; ++i) {
                /*
                 * only copy pointers
                 */
                Vector T1(_N, false, _T.column(i - 1));
                Vector T2(_N, false, _T.column(i - 2));

                dsbmv(2.0, _x, T1, 0., _T.column(i));
                daxpy(-1.0, T2, _T.column(i));
            }
        }
    }

    void getTransferFromParentChebyshevToChildrenChebyshev(index_t _nChebyshev, Vector &_chebyNode, Matrix &_tNode,
                                                           Matrix &_transfer) {
        Vector childChebyNode(2 * _nChebyshev);
        Vector T1(_nChebyshev);
        setValue(T1, -0.5);
        daxpy(0.5, _chebyNode, T1);
        memcpy(childChebyNode.data(), T1.data(), _nChebyshev * sizeof(scalar_t));

        Vector T2(_nChebyshev);
        setValue(T2, 0.5);
        daxpy(0.5, _chebyNode, T2);
        memcpy(childChebyNode.data() + _nChebyshev, T2.data(), _nChebyshev * sizeof(scalar_t));

        getStandardChebyPoly(_nChebyshev, 2 * _nChebyshev, childChebyNode, _transfer);

        Matrix T3(2 * _nChebyshev, _nChebyshev);
        setValue(T3, -1.0);

        dgemm_t(2.0, _transfer, _tNode, 1.0, T3);
        dscal(1.0 / _nChebyshev, T3);
        _transfer = T3;
    }


    void getTransfer(index_t _nChebyshev, Vector &_chebyNode, Matrix &_tNode, Matrix *R) {
        Matrix S(2 * _nChebyshev, _nChebyshev);
        getTransferFromParentChebyshevToChildrenChebyshev(_nChebyshev, _chebyNode, _tNode, S);

        Matrix Transfer[2];
        Transfer[0].resize(_nChebyshev, _nChebyshev);
        Transfer[1].resize(_nChebyshev, _nChebyshev);

        setBlock(Transfer[0], S, 0, 0, _nChebyshev, _nChebyshev);
        setBlock(Transfer[1], S, _nChebyshev, 0, _nChebyshev, _nChebyshev);

        index_t _rank = _nChebyshev * _nChebyshev;
        for (index_t i = 0; i < 4; ++i) {
            R[i].resize(_rank, _rank);
        }

        // follow bit representaion.
        for (index_t i = 0; i < _nChebyshev; ++i) {
            for (index_t j = 0; j < _nChebyshev; ++j) {
                for (index_t k = 0; k < _nChebyshev; ++k) {
                    for (index_t l = 0; l < _nChebyshev; ++l) {
                        for (index_t id = 0; id < 4; ++id) {
                            index_t bit[2];
                            bit[0] = (id >> 0) & 1;
                            bit[1] = (id >> 1) & 1;
                            R[id](i * _nChebyshev + j, k * _nChebyshev + l) =
                                    Transfer[bit[1]](i, k) * Transfer[bit[0]](j, l);
                        }
                    }
                }
            }
        }
    }


    void getScaledChebyNode(index_t _nChebyNode, Vector &_chebyNode, point &center, point &radius,
                            vector<point> &_scaledCnode) {
        for (index_t i = 0; i < _nChebyNode; ++i) {
            _scaledCnode.emplace_back(center.x + radius.x * _chebyNode(i),
                                         center.y + radius.y * _chebyNode(i));
        }
    }

    void getCharge(index_t rootId) {
        node &n = t.dict[rootId];
        if (n.chargeComputed) {
            return;
        } else {
            n.chargeComputed = true;
            n.charge.resize(n.nSource);
            for (index_t k = 0; k < n.nSource; ++k) {
                n.charge(k) = chargeTree(n.sourceIndex[k]);
            }
        }
    }

    void
    getTransferParentToChildren(index_t _nChebyNode, vector<point> &_tree, vector<index_t> &_index, point &_center,
                                point &_radius,
                                Vector &_chebyNode, Matrix &_tNode, Matrix &R) {

        auto N = (index_t) _index.size();
        Vector standlocation[DIM];
        standlocation[0].resize(N);
        standlocation[1].resize(N);
        for (index_t i = 0; i < N; ++i) {
            standlocation[0](i) = (_tree[_index[i]].x - _center.x) / _radius.x;
            standlocation[1](i) = (_tree[_index[i]].y - _center.y) / _radius.y;
        }

        Matrix Transfer[DIM];
        for (index_t k = 0; k < DIM; ++k) {
            Transfer[k].resize(N, _nChebyNode);
            getStandardChebyPoly(_nChebyNode, N, standlocation[k], Transfer[k]);

            Matrix T3(N, _nChebyNode);
            setValue(T3, -1.0);
            dgemm_t(2.0, Transfer[k], _tNode, 1.0, T3);
            dscal(1.0 / _nChebyNode, T3);
            Transfer[k] = T3;

        }

        index_t _rank = _nChebyNode * _nChebyNode;
        R.resize(N, _rank);
        for (index_t k = 0; k < N; ++k) {
            for (index_t i = 0; i < _nChebyNode; ++i) {
                for (index_t j = 0; j < _nChebyNode; ++j) {
                    R(k, j * _nChebyNode + i) = Transfer[0](k, i) * Transfer[1](k, j);
                }
            }
        }
    }

    void upPass(index_t rootId = 0) {
        node &n = t.dict[rootId];
        n.scaledCnode.clear();
        n.nodeCharge.resize(rank);
        n.nodePotential.resize(rank);
        getScaledChebyNode(nChebyshev, chebyNode, n.center, n.radius, n.scaledCnode);

        if (n.isLeaf) {
            // lazy
            getCharge(rootId);
            getTransferParentToChildren(nChebyshev, t.sourceTree, n.sourceIndex,
                    n.center, n.radius, chebyNode, tNode, n.R);
            /*
             * in case leaf node has no points, which causes DGEMV error.
             */
            if (n.R.row() != 0) dgemv_t(1.0, n.R, n.charge, 1.0, n.nodeCharge);
        } else {
            for (index_t i = 0; i < 4; ++i) {
#ifdef RUN_OMP
#pragma omp task shared(n) firstprivate(i)
#endif
                upPass(n.child[i]);
            }
#ifdef RUN_OMP
#pragma omp taskwait
#endif
            for (index_t i = 0; i < 4; ++i) {
                if (!t.dict[n.child[i]].isEmpty) {
                    dgemv_t(1.0, R[i], t.dict[n.child[i]].nodeCharge, 1.0, n.nodeCharge);
                }
            }
        }
    }

};

#endif //TREECODE_RTE_H