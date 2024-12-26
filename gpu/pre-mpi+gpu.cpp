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



struct Point {
    double x, y;

    Point(double x_val = 0.0, double y_val = 0.0) : x(x_val), y(y_val) {}
};

Point P1 = Point(0.0, 0.0);
Point P2 = Point(0.0, 3.0);
Point P3 = Point(3.0, 0.0);
Point P4 = Point(3.0, 3.0);




double stop = (9e-8);


typedef struct {
    int left, right, top, bottom;
} Neighbors;


typedef struct {
    double* left, *right, *top, *bottom;
} Boundary;


double* gen_grid(double start, double delta, int size) {
    double* res = (double*)malloc(sizeof(double) * size);
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

double gen_a_ij(int i, int j, Grid* grid, double h1, double h2) {

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
    if (y_2 < upb) upb = y_2;
    // upb = std::min(upb, y_2);
    double downb = y_1;

    double length = upb - downb;
    return length / h1 + (1.0 - length / h1) / grid->epsilon;
}

double gen_b_ij(int i, int j, Grid* grid, double h1, double h2) {
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
    if (x_2 < inter_x) inter_x = x_2;
    // inter_x = std::min(inter_x, x_2);

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

double gen_f_ij(int i, int j, Grid* grid, double h1, double h2) {
    
    
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






int main (int argc, char** argv) {
    MPI_Comm comm;
    int rank, size;
    int coord[2];

    int period[2], reorder;

    int M = 160;
    int N = 180;

    double h1 = (P4.y - P1.y) / N;
    double h2 = (P4.x - P1.x) / M;


    double scals[2] = {0, 0};
    double pretau[2] = {0, 0};

    MPI_Init(&argc, &argv);

    
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);


    int* dim = new int[2];
    dim[0] = 0;
    dim[1] = 0;

    MPI_Dims_create(size, 2, dim);


    period[0] = 0;
    period[1] = 0;
    reorder = 0;

    

    MPI_Cart_create(MPI_COMM_WORLD, 2, dim, period, reorder, &comm);

    MPI_Cart_coords(comm, rank, 2, coord);




    
    Neighbors neighbors = Neighbors({-1, -1, -1, -1});


    MPI_Cart_shift(comm, 0, 1, &neighbors.left, &neighbors.right);
    MPI_Cart_shift(comm, 1, 1, &neighbors.bottom, &neighbors.top);

    
    int dy = (N + 1) / dim[1];     // размер сетки одного процесса dim_1  - y
    int dx = (M + 1) / dim[0];     // размер сетки однного процесса dim_0 - x

    int offy = (N + 1) % dy;
    int offx = (M + 1) % dx;


    if (coord[1] < offy) {
        offy = coord[1];
    }
    if (coord[0] < offx) {
        offx = coord[0];
    }   
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



    Grid grid = Grid(P, h2, h1, cdx, cdy);




    double* operator_a = (double*)malloc(sizeof(double)*cdx*cdy);
    double* operator_b = (double*)malloc(sizeof(double)*cdx*cdy);
    double* operator_f = (double*)malloc(sizeof(double)*cdx*cdy);
    for (int i = 0; i < cdy; ++i) {
        for (int j = 0; j < cdx; ++j) {
            operator_a[i*cdx + j] = gen_a_ij(i, j, &grid, h1, h2);
            operator_b[i*cdx + j] = gen_b_ij(i, j, &grid, h1, h2);
            operator_f[i*cdx + j] = gen_f_ij(i, j, &grid, h1, h2);
        }
    }


    double* w = new double[cdx*cdy]();
    double* r = new double[cdx*cdy]();
    double* Ar = new double[cdx*cdy]();
    
    double* send_row = new double[cdy]();
    double* send_col = new double[cdx]();


    // Boundary ba, bb, bw, br;

    double* ba_top = new double[cdx]();
    double* ba_bottom = new double[cdx]();
    double* ba_left = new double[cdy]();
    double* ba_right = new double[cdy]();

    double* bb_top = new double[cdx]();
    double* bb_bottom = new double[cdx]();
    double* bb_left = new double[cdy]();
    double* bb_right = new double[cdy]();

    double* bw_top = new double[cdx]();
    double* bw_bottom = new double[cdx]();
    double* bw_left = new double[cdy]();
    double* bw_right = new double[cdy]();

    double* br_top = new double[cdx]();
    double* br_bottom = new double[cdx]();
    double* br_left = new double[cdy]();
    double* br_right = new double[cdy]();

    
    int it = 0;
    double error=999;

    int request_count = 0;

   


    MPI_Request requests0[8];
        request_count = 0;
    #pragma acc data copy(send_col[0:cdx], send_row[0:cdy],  w[0:cdx * cdy], r[0:cdx * cdy], Ar[0:cdx * cdy], operator_a[0:cdx * cdy], operator_b[0:cdx * cdy], operator_f[0:cdx * cdy], ba_top[0:cdx], ba_bottom[0:cdx], ba_left[0:cdy], ba_right[0:cdy], bb_top[0:cdx], bb_bottom[0:cdx], bb_left[0:cdy], bb_right[0:cdy], bw_top[0:cdx], bw_bottom[0:cdx], bw_left[0:cdy], bw_right[0:cdy], br_top[0:cdx], br_bottom[0:cdx], br_left[0:cdy], br_right[0:cdy])
    {
        if (neighbors.left >= 0) {
            #pragma acc parallel loop
            for (int i = 0; i < cdy; ++i) {
                send_row[i] = operator_a[i * cdx];
            }
            MPI_Isend(&(send_row), cdy, MPI_DOUBLE, neighbors.left, 0, comm, &requests0[request_count++]);
            MPI_Irecv(&(ba_left), cdy, MPI_DOUBLE, neighbors.left, 1, comm, &requests0[request_count++]);

            
        }

        if (neighbors.right >= 0) {
            #pragma acc parallel loop
            for (int i = 0; i < cdy; ++i) {
                
                send_row[i] = operator_a[i * (cdx) + (cdx - 1)];
            }
            MPI_Isend(&send_row, cdy, MPI_DOUBLE, neighbors.right, 1, comm, &requests0[request_count++]);
            MPI_Irecv(&(ba_right), cdy, MPI_DOUBLE, neighbors.right, 0, comm, &requests0[request_count++]);
            
        }

        if (neighbors.top >= 0) {
            #pragma acc parallel loop
            for (int j = 0; j < cdx; ++j) {
                send_col[j] = operator_a[(cdy-1) * cdx + j];
            }
            MPI_Isend(&send_col, M, MPI_DOUBLE, neighbors.top, 2, comm, &requests0[request_count++]);
            MPI_Irecv(&(ba_top), M, MPI_DOUBLE, neighbors.top, 3, comm, &requests0[request_count++]);
        }

        if (neighbors.bottom >= 0) {
            #pragma acc parallel loop
            for (int j = 0; j < cdx; ++j) {
                send_col[j] = operator_a[j];
            }

            MPI_Isend(&send_col, cdx, MPI_DOUBLE, neighbors.bottom, 3, comm, &requests0[request_count++]);
            MPI_Irecv(&(ba_bottom), cdx, MPI_DOUBLE, neighbors.bottom, 2, comm, &requests0[request_count++]);

        }

        MPI_Waitall(request_count, requests0, MPI_STATUSES_IGNORE);
    }

    MPI_Request requests01[8];
        request_count = 0;
    #pragma acc data copy(send_col[0:cdx], send_row[0:cdy],  w[0:cdx * cdy], r[0:cdx * cdy], Ar[0:cdx * cdy], operator_a[0:cdx * cdy], operator_b[0:cdx * cdy], operator_f[0:cdx * cdy], ba_top[0:cdx], ba_bottom[0:cdx], ba_left[0:cdy], ba_right[0:cdy], bb_top[0:cdx], bb_bottom[0:cdx], bb_left[0:cdy], bb_right[0:cdy], bw_top[0:cdx], bw_bottom[0:cdx], bw_left[0:cdy], bw_right[0:cdy], br_top[0:cdx], br_bottom[0:cdx], br_left[0:cdy], br_right[0:cdy])
    {
        if (neighbors.left >= 0) {
            #pragma acc parallel loop
            for (int i = 0; i < cdy; ++i) {
                send_row[i] = operator_b[i * cdx];
            }
            MPI_Isend(&(send_row), cdy, MPI_DOUBLE, neighbors.left, 0, comm, &requests01[request_count++]);
            MPI_Irecv(&(bb_left), cdy, MPI_DOUBLE, neighbors.left, 1, comm, &requests01[request_count++]);

            
        }

        if (neighbors.right >= 0) {
            #pragma acc parallel loop
            for (int i = 0; i < cdy; ++i) {
                
                send_row[i] = operator_b[i * (cdx) + (cdx - 1)];
            }
            MPI_Isend(&send_row, cdy, MPI_DOUBLE, neighbors.right, 1, comm, &requests01[request_count++]);
            MPI_Irecv(&(bb_right), cdy, MPI_DOUBLE, neighbors.right, 0, comm, &requests01[request_count++]);
            
        }

        if (neighbors.top >= 0) {
            #pragma acc parallel loop
            for (int j = 0; j < cdx; ++j) {
                send_col[j] = operator_b[(cdy-1) * cdx + j];
            }
            MPI_Isend(&send_col, M, MPI_DOUBLE, neighbors.top, 2, comm, &requests01[request_count++]);
            MPI_Irecv(&(bb_top), M, MPI_DOUBLE, neighbors.top, 3, comm, &requests01[request_count++]);
        }

        if (neighbors.bottom >= 0) {
            #pragma acc parallel loop
            for (int j = 0; j < cdx; ++j) {
                send_col[j] = operator_b[j];
            }

            MPI_Isend(&send_col, cdx, MPI_DOUBLE, neighbors.bottom, 3, comm, &requests01[request_count++]);
            MPI_Irecv(&(bb_bottom), cdx, MPI_DOUBLE, neighbors.bottom, 2, comm, &requests01[request_count++]);

        }

        MPI_Waitall(request_count, requests01, MPI_STATUSES_IGNORE);
    }



    double scal1 = 0;
    double scal2 = 0;

    double start = MPI_Wtime();
    
    #pragma acc data copy(send_col[0:cdx], send_row[0:cdy],  w[0:cdx * cdy], r[0:cdx * cdy], Ar[0:cdx * cdy], operator_a[0:cdx * cdy], operator_b[0:cdx * cdy], operator_f[0:cdx * cdy], ba_top[0:cdx], ba_bottom[0:cdx], ba_left[0:cdy], ba_right[0:cdy], bb_top[0:cdx], bb_bottom[0:cdx], bb_left[0:cdy], bb_right[0:cdy], bw_top[0:cdx], bw_bottom[0:cdx], bw_left[0:cdy], bw_right[0:cdy], br_top[0:cdx], br_bottom[0:cdx], br_left[0:cdy], br_right[0:cdy])
    while(error > stop){

        MPI_Request requests[8];
        request_count = 0;
        if (neighbors.left >= 0) {
            #pragma acc parallel loop
            for (int i = 0; i < cdy; ++i) {
                send_row[i] = w[i * cdx];
            }
            MPI_Isend(&(send_row), cdy, MPI_DOUBLE, neighbors.left, 0, comm, &requests[request_count++]);
            MPI_Irecv(&(bw_left), cdy, MPI_DOUBLE, neighbors.left, 1, comm, &requests[request_count++]);

            
        }

        if (neighbors.right >= 0) {
            #pragma acc parallel loop
            for (int i = 0; i < cdy; ++i) {
                
                send_row[i] = w[i * (cdx) + (cdx - 1)];
            }
            MPI_Isend(&send_row, cdy, MPI_DOUBLE, neighbors.right, 1, comm, &requests[request_count++]);
            MPI_Irecv(&(bw_right), cdy, MPI_DOUBLE, neighbors.right, 0, comm, &requests[request_count++]);
            
        }

        if (neighbors.top >= 0) {
            #pragma acc parallel loop
            for (int j = 0; j < cdx; ++j) {
                send_col[j] = w[(cdy-1) * cdx + j];
            }
            MPI_Isend(&send_col, M, MPI_DOUBLE, neighbors.top, 2, comm, &requests[request_count++]);
            MPI_Irecv(&(bw_top), M, MPI_DOUBLE, neighbors.top, 3, comm, &requests[request_count++]);
        }

        if (neighbors.bottom >= 0) {
            #pragma acc parallel loop
            for (int j = 0; j < cdx; ++j) {
                send_col[j] = w[j];
            }

            MPI_Isend(&send_col, cdx, MPI_DOUBLE, neighbors.bottom, 3, comm, &requests[request_count++]);
            MPI_Irecv(&(bw_bottom), cdx, MPI_DOUBLE, neighbors.bottom, 2, comm, &requests[request_count++]);

        }

        MPI_Waitall(request_count, requests, MPI_STATUSES_IGNORE);


        
        #pragma acc parallel loop collapse(2)
        for (int i = 0; i < cdy; ++i) {
            for (int j = 0; j < cdx; ++j) {
                if (i + offy == 0 || i + offy == N || j + offx == 0 || j + offx == M) {
                    r[i*cdx + j] = 0;
                    continue;
                }
                double z[6];
                double x[6];
                
                if (i - 1 < 0) {
                    x[5] = bw_bottom[j];
                } else {
                    x[5] = w[(i-1)*cdx + j];
                }

                if (i + 1>= cdy) {
                    x[1] = bw_top[j];
                    x[0] = bb_top[j];

                } else {
                    x[1] = w[(i+1)*cdx + j];
                    x[0] = operator_b[(i+1)*cdx + j];
                }
                
                if (j - 1 < 0) {
                    z[5] = bw_left[i];
                } else {
                    z[5] = w[i*cdx + (j-1)];
                }
                if (j + 1 >= cdx) {
                    z[0] = ba_right[i];
                    z[1] = bw_right[i];
                } else {
                    z[0] = operator_a[i*cdx + (j+1)];
                    z[1] = w[i*cdx + (j+1)];
                }

                z[2] = w[i*cdx + j];
                z[3] = operator_a[i*cdx + j];
                z[4] = w[i*cdx + j];
                x[2] = w[i*cdx + j];
                x[3] = operator_b[i*cdx + j];
                x[4] = w[i*cdx + j];
                
                r[i*cdx + j] = 
                    (-1.0 / h2) * 
                    (
                        z[0] * 
                        (z[1] - z[2]) / h2 - 
                        z[3] * 
                        (z[4]  - z[5]) / h2
                    ) +
                    (-1.0 / h1) * 
                    (
                        x[0] * 
                        (x[1] - x[2]) / h1 - 
                        x[3] * 
                        (x[4] - x[5]) / h1
                    );
                
            }   
        }
        #pragma acc parallel loop collapse(2)    
        for (int i = 0; i < cdy; ++i) {
            for (int j = 0; j < cdx; ++j) {
                if ((offx + j != 0) && (offx + j != (M)) && (offy + i != 0) && (offy + i != (N))){
                    r[i*cdx + j] -= operator_f[i*cdx + j];
                }else {
                    r[i*cdx + j] = 0;
                }
            }
        }

        

        #pragma acc update self(br_left[0:cdy])
        #pragma acc update self(br_right[0:cdy])
        #pragma acc update self(br_top[0:cdx])
        #pragma acc update self(br_bottom[0:cdx])
        
        MPI_Request requests1[8];
        request_count = 0;

        if (neighbors.left >= 0) {
            #pragma acc parallel loop
            for (int i = 0; i < cdy; ++i) {
                send_row[i] = r[i * cdx];
            }
            MPI_Isend(&(send_row), cdy, MPI_DOUBLE, neighbors.left, 0, comm, &requests1[request_count++]);
            MPI_Irecv(&(br_left), cdy, MPI_DOUBLE, neighbors.left, 1, comm, &requests1[request_count++]);

            
        }

        if (neighbors.right >= 0) {
            #pragma acc parallel loop
            for (int i = 0; i < cdy; ++i) {
                
                send_row[i] = r[i * (cdx) + (cdx - 1)];
            }
            MPI_Isend(&send_row, cdy, MPI_DOUBLE, neighbors.right, 1, comm, &requests1[request_count++]);
            MPI_Irecv(&(br_right), cdy, MPI_DOUBLE, neighbors.right, 0, comm, &requests1[request_count++]);
            
        }

        if (neighbors.top >= 0) {
            #pragma acc parallel loop
            for (int j = 0; j < cdx; ++j) {
                send_col[j] = r[(cdy-1) * cdx + j];
            }
            MPI_Isend(&send_col, M, MPI_DOUBLE, neighbors.top, 2, comm, &requests1[request_count++]);
            MPI_Irecv(&(br_top), M, MPI_DOUBLE, neighbors.top, 3, comm, &requests1[request_count++]);
        }

        if (neighbors.bottom >= 0) {
            #pragma acc parallel loop
            for (int j = 0; j < cdx; ++j) {
                send_col[j] = r[j];
            }

            MPI_Isend(&send_col, cdx, MPI_DOUBLE, neighbors.bottom, 3, comm, &requests1[request_count++]);
            MPI_Irecv(&(br_bottom), cdx, MPI_DOUBLE, neighbors.bottom, 2, comm, &requests1[request_count++]);

        }


        MPI_Waitall(request_count, requests1, MPI_STATUSES_IGNORE);



        #pragma acc update device(br_left[0:cdy])
        #pragma acc update device(br_right[0:cdy])
        #pragma acc update device(br_top[0:cdx])
        #pragma acc update device(br_bottom[0:cdx])

        #pragma acc parallel loop collapse(2)
        for (int i = 0; i < cdy; ++i) {
            for (int j = 0; j < cdx; ++j) {
                if (i + offy == 0 || i + offy == N || j + offx == 0 || j + offx == M) {
                    Ar[i*cdx + j] = 0;
                    continue;
                }

                double z[6];
                double x[6];
                
                if (i - 1 < 0) {
                    x[5] = br_bottom[j];
                } else {
                    x[5] = r[(i-1)*cdx + j];
                }

                if (i + 1>= cdy) {
                    x[1] = br_top[j];
                    x[0] = bb_top[j];

                } else {
                    x[1] = r[(i+1)*cdx + j];
                    x[0] = operator_b[(i+1)*cdx + j];
                }
                
                if (j - 1 < 0) {
                    z[5] = br_left[i];
                } else {
                    z[5] = r[i*cdx + (j-1)];
                }
                if (j + 1 >= cdx) {
                    z[0] = ba_right[i];
                    z[1] = br_right[i];
                } else {
                    z[0] = operator_a[i*cdx + (j+1)];
                    z[1] = r[i*cdx + (j+1)];
                }
                z[2] = r[i*cdx + j];
                z[3] = operator_a[i*cdx + j];
                z[4] = r[i*cdx + j];
                x[2] = r[i*cdx + j];
                x[3] = operator_b[i*cdx + j];
                x[4] = r[i*cdx + j];
                
                Ar[i*cdx + j] = 
                    (-1.0 / h2) * 
                    (
                        z[0] * 
                        (z[1] - z[2]) / h2 - 
                        z[3] * 
                        (z[4]  - z[5]) / h2
                    ) +
                    (-1.0 / h1) * 
                    (
                        x[0] * 
                        (x[1] - x[2]) / h1 - 
                        x[3] * 
                        (x[4] - x[5]) / h1
                    );
            }   
        }


        scal1 = 0;
        #pragma acc parallel loop collapse(2)  reduction(+:scal1)
        for (int i = 0; i < cdy; ++i) {
            for (int j = 0; j < cdx; ++j) {
                scal1 += h2 * h1 * r[i*cdx + j] * r[i*cdx + j];
            }
        }
        

        scal2 = 0;
        #pragma acc parallel loop collapse(2)  reduction(+:scal2)
        for (int i = 0; i < cdy; ++i) {
            for (int j = 0; j < cdx; ++j) {
                scal2 += h2 * h1 * Ar[i*cdx + j] * r[i*cdx + j];
            }
        } 

        scals[0] = scal1;
        scals[1] = scal2;
	    MPI_Allreduce(scals, pretau, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        double tau = pretau[0] / pretau[1];
        

        #pragma acc parallel loop collapse(2)
        for (int i = 0; i < cdy; ++i) {
            for (int j = 0; j < cdx; ++j) {
                w[i*cdx + j] -= tau * r[i*cdx + j];
            }
        }   
        
        error = sqrt(tau * tau * pretau[0]);
        it += 1;
    }

    if (rank == 0) {
        printf("%d %f\n", it, MPI_Wtime() - start);
    }

    MPI_Finalize();

    return 0;
}
