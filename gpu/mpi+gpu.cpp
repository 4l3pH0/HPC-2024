#include <iostream>
#include <fstream>
#include <stdio.h>
#include <cmath>
#include <math.h>
#include <mpi.h>
#include <omp.h>
#include <iomanip>
#include <vector>
#include <chrono>

// var 4


const int M = 160;
const int N = 180;

struct Point {
    double x, y;

    Point(double x_val = 0.0, double y_val = 0.0) : x(x_val), y(y_val) {}
};

Point P1 = Point(0.0, 0.0);
Point P2 = Point(0.0, 3.0);
Point P3 = Point(3.0, 0.0);
Point P4 = Point(3.0, 3.0);

double h1 = (P4.y - P1.y) / N;
double h2 = (P4.x - P1.x) / M;



double stop = (9e-8);


typedef struct {
    int left, right, top, bottom;
} Neighbors;


typedef struct {
    double* left, right, top, bottom;
} Boundary;



double* gen_grid(double start, double delta, int size) {
    double* res(size, 0);
    #pragma omp parallel for
    for (int i = 0; i < size; ++i) {
        res[i] = start + i * delta;
    }
    return res; 
}

struct Grid {
    double* dx;
    double* dy;
    double epsilon;

    Grid(Point start, double h2, double h1, int M, int N) 
        : dx(gen_grid(start.x, h2, M)), 
          dy(gen_grid(start.y, h1, N)) 
    {
        epsilon = std::max(h2,h1);
        epsilon *= epsilon;
    }
};

double cross_product(const Point& p0, const Point& p1, const Point& p2) {
    return (p1.x - p0.x) * (p2.y - p0.y) - (p1.y - p0.y) * (p2.x - p0.x);
}

struct Operator {
    double* operator_a;
    double* operator_b;
    double* operator_f;
    Grid* grid;

    Operator(const double*& operator_a,
             const double*& operator_b,
             const double*& operator_f,
             Grid* grid)
        : operator_a(operator_a), operator_b(operator_b), operator_f(operator_f), grid(grid) {}
};

bool is_point_in_trapezoid(const Point& p) {
    Point A = Point(0, 0);
    Point B = Point(3, 0);
    Point C = Point(2, 3);
    Point D = Point(0, 3);

    bool on_left_of_AB = cross_product(A, B, p) >= 0;
    bool on_left_of_BC = cross_product(B, C, p) >= 0;
    bool on_left_of_CD = cross_product(C, D, p) >= 0;
    bool on_left_of_DA = cross_product(D, A, p) >= 0;

    return on_left_of_AB && on_left_of_BC && on_left_of_CD && on_left_of_DA;
}

double gen_a_ij(int i, int j, Grid* grid) {

    double x = grid->dx[j] - 0.5 * h2;
    double y_1 = grid->dy[i] - 0.5 * h1;
    double y_2 = grid->dy[i] + 0.5 * h1;
    Point p1(x, y_1);
    Point p2(x, y_2);
    if (is_point_in_trapezoid(p1) && is_point_in_trapezoid(p2)) {
        return 1.0;
    }
    if (!is_point_in_trapezoid(p1) && !is_point_in_trapezoid(p2)) {
        return 1 / grid->epsilon;
    }

    double upb = (x <= 2.0) ? 3.0 : -3 * x + 9.0;
    upb = std::min(upb, y_2);
    double downb = y_1;

    double length = upb - downb;
    return length / h1 + (1.0 - length / h1) / grid->epsilon;
}

double gen_b_ij(int i, int j, Grid* grid) {
    double y = grid->dy[i] - 0.5 * h1;

    double x_1 = grid->dx[j] - 0.5 * h2;
    double x_2 = grid->dx[j] + 0.5 * h2;
    Point p1(x_1, y);
    Point p2(x_2, y);
    if (is_point_in_trapezoid(p1) && is_point_in_trapezoid(p2)) {
        return 1.0;
    }
    if (!is_point_in_trapezoid(p1) && !is_point_in_trapezoid(p2)) {
        return 1 / grid->epsilon;
    }

    double inter_x = 3.0 - y / 3.0;
    inter_x = std::min(inter_x, x_2);

    double upb = inter_x;
    double downb = x_1;
    double length = upb - downb;
    return length / h2 + (1.0 - length / h2) / grid->epsilon;
}

Point intersect(Point& p1, Point& p2, Point& p3, Point& p4) {
    double A1 = p2.y - p1.y;
    double B1 = p1.x - p2.x;
    double C1 = A1 * p1.x + B1 * p1.y;

    double A2 = p4.y - p3.y;
    double B2 = p3.x - p4.x;
    double C2 = A2 * p3.x + B2 * p3.y;
    double det = A1 * B2 - A2 * B1;
    double x = (B2 * C1 - B1 * C2) / det;
    double y = (A1 * C2 - A2 * C1) / det;
    return Point{x, y};
          
} 

double gen_f_ij(int i, int j, Grid* grid) {
    
    
    double y_1 = grid->dy[i] - 0.5 * h1;
    double y_2 = grid->dy[i] + 0.5 * h1;    

    double x_1 = grid->dx[j] - 0.5 * h2;
    double x_2 = grid->dx[j] + 0.5 * h2;

    if (x_1 > 3) return 0;
    if (x_2 < 0) return 0;
    if (y_1 > 3) return 0;
    if (y_2 < 0) return 0;

    if (x_2 > 3) x_2 = 3;
    if (y_2 > 3) y_2 = 3;
    if (x_1 < 0) x_1 = 0;
    if (y_1 < 0) y_1 = 0;

    Point A_left_down(x_1, y_1);
    Point A_left_up(x_1, y_2);
    Point A_right_up(x_2, y_2);
    Point A_right_down(x_2, y_1);
    Point p1(3, 0);
    Point p2(2, 3);

    
    int cnt = int(is_point_in_trapezoid(A_left_down)) + int(is_point_in_trapezoid(A_left_up)) + int(is_point_in_trapezoid(A_right_up)) + int(is_point_in_trapezoid(A_right_down));

    if (cnt == 4) {
        return (x_2 - x_1) * (y_2 - y_1) / (h2 * h1);
    }
    if (cnt == 3) {
        Point inter1 = intersect(A_left_up, A_right_up, p1, p2);
        Point inter2 = intersect(A_right_down, A_right_up, p1, p2);
        double a = x_2 - inter1.x;
        double b = y_2 - inter2.y;
        return ((x_2 - x_1)*(y_2 - y_1) - a*b/2)/ (h2 * h1);
    }
    if (cnt == 2) {
        if (is_point_in_trapezoid(A_left_up) ){
            Point inter1 = intersect(A_left_up, A_right_up, p1, p2);
            Point inter2 = intersect(A_left_down, A_right_down, p1, p2);    
	    double h = y_2-y_1;
            double s1 = inter1.x - x_1;
            double s2 = inter2.x - x_1;
            return (0.5 * h * (s1 + s2) / 2) / (h2 * h1);
        }
        if (is_point_in_trapezoid(A_right_down)){
            double h = x_2 - x_1;
            Point inter1 = intersect(A_left_down, A_left_up, p1, p2);
	    Point inter2 = intersect(A_right_down, A_right_up, p1, p2);
	    double s1 = inter1.y - y_1;
            double s2 = inter2.y - y_1;
            return (0.5 * h * (s1 + s2) / 2) / (h2 * h1);
        }
    }
    if (cnt == 1) {
        Point inter1 = intersect(A_left_down, A_left_up, p1, p2);
	Point inter2 = intersect(A_left_down, A_right_down, p1, p2);
        double a = inter1.y - y_1;
        double b = inter2.y - x_1;
        return 0.5 *a*b/  (h2 * h1);
    }
    if (cnt == 0) {
        return 0;
    }

    return 0;
}

Operator* do_op(Grid* grid, int dx, int dy) {
    double* operator_a(dx*dy, 0);
    double* operator_b(dx*dy, 0);
    double* operator_f(dx*dy, 0);
    #pragma omp parallel for collapse(2)
    for (int i = 0; i < dy; ++i) {
        for (int j = 0; j < dx; ++j) {
            operator_a[i*dx + j] = gen_a_ij(i, j, grid);
            operator_b[i*dx + j] = gen_b_ij(i, j, grid);
            operator_f[i*dx + j] = gen_f_ij(i, j, grid);
        }
    }
    return new Operator(operator_a, operator_b, operator_f, grid);

}

void exchange_boundaries(double*& grid, const Neighbors& neighbors, MPI_Comm& comm, Boundary& bond, double*& send_col, double*& send_row,  int& M, int& N) {

    MPI_Request requests[8];
    int request_count = 0;
    // MPI_Status status;

    // double* left, right, top, bottom;

    // Exchange with left neighbor
    if (neighbors.left >= 0) {
        #pragma omp parallel for
        for (int i = 0; i < N; ++i) {
            send_row[i] = grid[i * M];
        }
        MPI_Isend(send_row.data(), N, MPI_DOUBLE, neighbors.left, 0, comm, &requests[request_count++]);
        MPI_Irecv(bond.left.data(), N, MPI_DOUBLE, neighbors.left, 1, comm, &requests[request_count++]);

        
    }

    // Exchange with right neighbor
    if (neighbors.right >= 0) {
        #pragma omp parallel for
        for (int i = 0; i < N; ++i) {
            
            send_row[i] = grid[i * (M) + (M - 1)];
        }
        MPI_Isend(send_row.data(), N, MPI_DOUBLE, neighbors.right, 1, comm, &requests[request_count++]);
        MPI_Irecv(bond.right.data(), N, MPI_DOUBLE, neighbors.right, 0, comm, &requests[request_count++]);
        
    }

    // Exchange with top neighbor
    if (neighbors.top >= 0) {
        #pragma omp parallel for
        for (int j = 0; j < M; ++j) {
            send_col[j] = grid[(N-1) * M + j];
        }
        MPI_Isend(send_col.data(), M, MPI_DOUBLE, neighbors.top, 2, comm, &requests[request_count++]);
        MPI_Irecv(bond.top.data(), M, MPI_DOUBLE, neighbors.top, 3, comm, &requests[request_count++]);
    }

    // Exchange with bottom neighbor
    if (neighbors.bottom >= 0) {
        #pragma omp parallel for
        for (int j = 0; j < M; ++j) {
            send_col[j] = grid[j];
        }

        MPI_Isend(send_col.data(), M, MPI_DOUBLE, neighbors.bottom, 3, comm, &requests[request_count++]);
        MPI_Irecv(bond.bottom.data(), M, MPI_DOUBLE, neighbors.bottom, 2, comm, &requests[request_count++]);

    }

    MPI_Waitall(request_count, requests, MPI_STATUSES_IGNORE);
}

double choose_wisely(double*& a, int i, int j, int dx, int dy, Boundary& b) {
    if (i < 0) {
        return b.bottom[j];
    }
    if (i >= dy) {
        return b.top[j];
    }
    if (j < 0) {
        return b.left[i];
    }
    if (j >= dx) {
        return b.right[i];
    }
    return a[i*dx + j];
}

void solve_iter(double*& a, double*& b, double*& w, double*& r, int& dx, int& dy, int& offx, int& offy, Boundary& ba, Boundary& bb, Boundary& bw) {
    #pragma omp parallel for collapse(2)
    for (int i = 0; i < dy; ++i) {
        for (int j = 0; j < dx; ++j) {
            if (i + offy == 0 || i + offy == N || j + offx == 0 || j + offx == M) {
                r[i*dx + j] = 0;
                continue;
            }

            r[i*dx + j] = 
                (-1.0 / h2) * 
                (
                    choose_wisely(a, i, j + 1, dx, dy, ba) * 
                    (choose_wisely(w,i,j + 1,dx,dy,bw) - choose_wisely(w,i,j,dx,dy,bw)) / h2 - 
                    choose_wisely(a, i, j, dx, dy, ba) * 
                    (choose_wisely(w,i,j,dx,dy,bw)  - choose_wisely(w,i,j-1,dx,dy,bw)) / h2
                ) +
                (-1.0 / h1) * 
                (
                    choose_wisely(b, i + 1, j, dx, dy, bb) * 
                    (choose_wisely(w,i+1,j,dx,dy,bw) - choose_wisely(w,i,j,dx,dy,bw)) / h1 - 
                    choose_wisely(b, i, j, dx, dy, bb) * 
                    (choose_wisely(w,i,j,dx,dy,bw) - choose_wisely(w,i-1,j,dx,dy,bw)) / h1
                );
                // (-1.0 / h2) * 
                // (
                //     choose_wisely(a, i+1, j, dx, dy, ba) * 
                //     (choose_wisely(w,i+1,j,dx,dy,bw) - choose_wisely(w,i,j,dx,dy,bw)) / h2 - 
                //     choose_wisely(a, i, j, dx, dy, ba) * 
                //     (choose_wisely(w,i,j,dx,dy,bw)  - choose_wisely(w,i-1,j,dx,dy,bw)) / h2
                // ) +
                // (-1.0 / h1) * 
                // (
                //     choose_wisely(b, i, j+1, dx, dy, bb) * 
                //     (choose_wisely(w,i,j+1,dx,dy,bw) - choose_wisely(w,i,j,dx,dy,bw)) / h1 - 
                //     choose_wisely(b, i, j, dx, dy, bb) * 
                //     (choose_wisely(w,i,j,dx,dy,bw) - choose_wisely(w,i,j-1,dx,dy,bw)) / h1
                // );
        }   
    }
}

double scal_prod(double*& a, double*& b, int& dx, int& dy) {
    double sum = 0;
    #pragma omp parallel for collapse(2) reduction(+:sum)
    for (int i = 0; i < dy; ++i) {
        for (int j = 0; j < dx; ++j) {
            sum += h2 * h1 * a[i*dx + j] * b[i*dx + j];
        }
    }
    return sum;
}

void w_iter(double*& w, double*& r, double tau, int& dx, int& dy) {
    #pragma omp parallel for collapse(2)
    for (int i = 0; i < dy; ++i) {
        for (int j = 0; j < dx; ++j) {
            w[i*dx + j] -= tau * r[i*dx + j];
        }
    }
}

void init_bond(Boundary& bond, int dx, int dy) {
    bond.left.reserve(dy);
    bond.right.reserve(dy);
    bond.top.reserve(dx);
    bond.bottom.reserve(dx);
}



int main (int argc, char** argv) {
    #ifdef _OPENMP
        printf ("OpenMP is supported!\n");
    #endif
    MPI_Comm comm;
    int rank, size;
    int coord[2], id;

    int period[2], reorder;

    int request_count = 0;

    double scals[2];
    double pretau[2];

    std::vector<MPI_Request> requests(8);

    MPI_Init(&argc, &argv);

    
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);


    std::vector<int> dim(2, 0);

    MPI_Dims_create(size, 2, dim.data());


    period[0] = 0;
    period[1] = 0;
    reorder = 0;

    

    MPI_Cart_create(MPI_COMM_WORLD, 2, dim.data(), period, reorder, &comm);

    MPI_Cart_coords(comm, rank, 2, coord);




    
    Neighbors neighbors = Neighbors({-1, -1, -1, -1});


    MPI_Cart_shift(comm, 0, 1, &neighbors.left, &neighbors.right);
    MPI_Cart_shift(comm, 1, 1, &neighbors.bottom, &neighbors.top);

    
    int dy = (N + 1) / dim[1];     // размер сетки одного процесса dim_1  - y
    int dx = (M + 1) / dim[0];     // размер сетки однного процесса dim_0 - x

    int offy = (N + 1) % dy;
    int offx = (M + 1) % dx;

    // printf("qweqewwq   %d %d\n%d %d\n", coord[0], rank / dim[1], coord[1], rank % dim[1]);

    if (coord[1] < offy) {
        offy = coord[1];
    }
    if (coord[0] < offx) {
        offx = coord[0];
    }   
    //printf("%d %d %d", rank, coord[0], coord[1]);
    offy += dy * (coord[1]);
    offx += dx * (coord[0]);

    int cdx = dx;
    if (coord[0] < (M + 1) % dx) {
        ++cdx;
    }

    int cdy = dy;
    if (coord[1] < (N + 1) % dy) {
        ++cdy;
    }

    double px = P1.x + offx * h2;          // P1.x + dx*C_x
    double py = P1.y + offy * h1;          // P1.y + dy*C_y

    



    Point P = Point({px, py});


    // printf("point inited\n");


    Grid grid = Grid(P, h2, h1, cdx, cdy);

    // printf("grid inited\n");
    
    Operator* A = do_op(&grid, cdx, cdy);

    // printf("operator inited\n");

    double* w(cdx*cdy + 1, 0);
    double* r(cdx*cdy + 1, 0);
    double* Ar(cdx*cdy + 1, 0);
    
    double* send_row(cdy, 0);
    double* send_col(cdx, 0);


    Boundary ba, bb, bw, br;
    init_bond(ba, cdx, cdy);
    init_bond(bb, cdx, cdy);
    init_bond(bw, cdx, cdy);
    init_bond(br, cdx, cdy);

    double st;
    double end;
    
    int it = 0;
    double error=999;

   // double tau = 0;

   

   exchange_boundaries(A->operator_a, neighbors, comm, ba, send_col, send_row, cdx, cdy);
   exchange_boundaries(A->operator_b, neighbors, comm, bb, send_col, send_row, cdx, cdy);

    double start = MPI_Wtime();
    #pragma acc data copy(w[0:cdx * cdy], r[0:cdx * cdy], Ar[0:cdx * cdy], operatorA_a[0:cdx * cdy], operatorA_b[0:cdx * cdy], operatorA_f[0:cdx * cdy], send_row[0: cdy], send_col[0: cdx]])
    while(error > stop){
        // exchange_boundaries(w, neighbors, comm, bw, send_col, send_row, cdx, cdy);

        MPI_Request requests[8];
        int request_count = 0;
        // MPI_Status status;

        double* left, right, top, bottom;

        // Exchange with left neighbor
        if (neighbors.left >= 0) {
            for (int i = 0; i < cdy; ++i) {
                send_row[i] = w[i * cdx];
            }
            MPI_Isend(send_row.data(), N, MPI_DOUBLE, neighbors.left, 0, comm, &requests[request_count++]);
            MPI_Irecv(bw.left.data(), N, MPI_DOUBLE, neighbors.left, 1, comm, &requests[request_count++]);

            
        }

        // Exchange with right neighbor
        if (neighbors.right >= 0) {
            for (int i = 0; i < cdy; ++i) {
                
                send_row[i] = w[i * (cdx) + (cdx - 1)];
            }
            MPI_Isend(send_row.data(), N, MPI_DOUBLE, neighbors.right, 1, comm, &requests[request_count++]);
            MPI_Irecv(bw.right.data(), N, MPI_DOUBLE, neighbors.right, 0, comm, &requests[request_count++]);
            
        }

        // Exchange with top neighbor
        if (neighbors.top >= 0) {
            for (int j = 0; j < cdx; ++j) {
                send_col[j] = w[(cdy-1) * cdx + j];
            }
            MPI_Isend(send_col.data(), M, MPI_DOUBLE, neighbors.top, 2, comm, &requests[request_count++]);
            MPI_Irecv(bw.top.data(), M, MPI_DOUBLE, neighbors.top, 3, comm, &requests[request_count++]);
        }

        // Exchange with bottom neighbor
        if (neighbors.bottom >= 0) {
            for (int j = 0; j < cdx; ++j) {
                send_col[j] = w[j];
            }

            MPI_Isend(send_col.data(), M, MPI_DOUBLE, neighbors.bottom, 3, comm, &requests[request_count++]);
            MPI_Irecv(bw.bottom.data(), M, MPI_DOUBLE, neighbors.bottom, 2, comm, &requests[request_count++]);

        }

        MPI_Waitall(request_count, requests, MPI_STATUSES_IGNORE);




        solve_iter(A->operator_a, A->operator_b, w, r, cdx, cdy, offx, offy, ba, bb, bw);

        #pragma omp parallel for collapse(2)    
        for (int i = 0; i < cdy; ++i) {
            for (int j = 0; j < cdx; ++j) {
                if ((offx + j != 0) && (offx + j != (M)) && (offy + i != 0) && (offy + i != (N))){
                    r[i*cdx + j] -= A->operator_f[i*cdx + j];
                }else {
                    r[i*cdx + j] = 0;
                }
            }
        }

        
        exchange_boundaries(r, neighbors, comm, br, send_col, send_row, cdx, cdy);


        solve_iter(A->operator_a, A->operator_b, r, Ar, cdx, cdy, offx, offy, ba, bb, br);

        double scal1 = scal_prod(r, r, cdx, cdy);
        double scal2 = scal_prod(Ar, r, cdx, cdy);    

        scals[0] = scal1;
	scals[1] = scal2;
	MPI_Allreduce(scals, pretau, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        double tau = pretau[0] / pretau[1];
        
        // double tau = scal1/scal2; 
        w_iter(w, r, tau, cdx, cdy);
        
        error = sqrt(tau * tau * pretau[0]);
        it += 1;

     //    printf("%d %d %f %f\n", rank, it, error, tau);
    }

    if (rank == 0) {
        printf("%d %f\n", it, MPI_Wtime() - start);
    }

    MPI_Finalize();

    return 0;
}
