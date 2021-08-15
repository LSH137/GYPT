#include<stdio.h>
#include<stdlib.h>
#include "TH1.h"
#include "Riostream.h"
#include <string.h>
#include <dirent.h>
#include <math.h>
#include <algorithm>

#define TRUE 1
#define FALSE 0
#define MAX_DATA_LENGTH 300
#define MAX_EXPERIENCE_NUM 50
#define CM_HIGHT 0.038 //m
#define BOTTLE_HIGHT 0.205 //m
#define INITAL_FRAME_NUMBER 5 // must be > 1
#define TIME_INTERVER_ratio 0.133
#define MIN_TIME_INTERVER 0.0001
#define TOUCH_FLOOR 20
#define MIN_X -0.25
#define MAX_X 0.75
#define MIN_Y -0.3
#define MAX_Y 0.3
#define g 9.87

class Vector
{
public:
    Vector()
    {
        x = 0.0;
        y = 0.0;
        z = 0.0;
        t = 0.0;
    }
    ~Vector() = default;

    double x;
    double y;
    double z;
    double t;

    void reset()
    {
        x = 0.0;
        y = 0.0;
        z = 0.0;
        t = 0.0;
    }

    void print()
    {
        printf("%lfi + %lfj + %lfk, t = %lf, size = %lf\n", x, y, z, t, sqrt(x*x + y*y + z*z));
    }

    double getSize()
    {
        return sqrt(x*x + y*y + z*z);
    }

    void divide(double d)
    {
        if(d != 0)
        {
            x = x / d;
            y = y / d;
        }
        else
            printf("ERROR!: class: Vector method: divide() - divide by 0");
    }

};

bool compare_by_y(Vector v1, Vector v2)
{
    return v1.y < v2.y;
}

class Bottle
{
public:
    Bottle()
    {
        v_cm_i.reset();
        r_cm_i.reset();

        for(int i = 0; i < MAX_DATA_LENGTH; i++)
        {
            path_cm[i].reset();
            path_top[i].reset();
            w[i].reset();
        }
        for(int i = 0; i < INITAL_FRAME_NUMBER + 5; i++)
        {
            path_bottom[i].reset();
        }
        index_path_bottom = 0;
        index_path_cm = 0;
        index_path_top = 0;
        index_w = 0;
    }
    ~Bottle()
    {
        //delete[] path_cm;
        //delete[] path_top;
        //delete[] w;
    }

    Vector v_cm_i; // inital velocitu of center of mass
    Vector r_cm_i; // inital position of center of mass
    Vector* path_cm = new Vector[MAX_DATA_LENGTH];
    Vector* path_top = new Vector[MAX_DATA_LENGTH];
    Vector* w = new Vector[MAX_DATA_LENGTH]; // angular velocity(rad/s) by time
    Vector path_bottom[INITAL_FRAME_NUMBER + 5];

    void setPath_cm(Vector p)
    {
        if(index_path_cm < MAX_DATA_LENGTH - 1)
        {
            path_cm[index_path_cm] = p;
            index_path_cm++;
        }
        else
            printf("ERROR: location: Bottle > setPath_cm() - overflow error || index_path_cm = %d\n", index_path_cm);

    }

    void setPath_top(Vector p)
    {
        if(index_path_top < MAX_DATA_LENGTH - 1)
        {
            path_top[index_path_top] = p;
            index_path_top++;
            //printf("---- index_path_top = %d ----\n", index_path_top);
        }
        else
            printf("ERROR: location: Bottle > setPath_top() - overflow error || index_path_top = %d\n", index_path_top);

    }

    void setPath_bottom(Vector p)
    {
        if(index_path_bottom < MAX_DATA_LENGTH - 1)
        {
            path_bottom[index_path_bottom] = p;
            index_path_bottom++;
        }
        else
            printf("ERROR: location: Bottle > setPath_bottom() - overflow error || index_path_bottom = %d\n", index_path_bottom);
    }

    void setW(Vector p)
    {
        if(index_w < MAX_DATA_LENGTH - 1)
        {
            w[index_w] = p;
            index_w++;
        }
        else
            printf("ERROR: location: Bottle > setW() - overflow error\n");

    }

    double getMinY()
    {
        Vector* path_temp_top = new Vector[MAX_DATA_LENGTH];
        double min_y = 0.0;

        for(int i = 0; i < MAX_DATA_LENGTH; i++)
            path_temp_top[i] = path_top[i];
        
        std::sort(path_temp_top, path_temp_top + MAX_DATA_LENGTH, compare_by_y);
        min_y = path_temp_top[0].y;

        delete[] path_temp_top;

        return min_y;
    }

    void resetBottle()
    {
        v_cm_i.reset();
        r_cm_i.reset();
        for(int i = 0; i < MAX_DATA_LENGTH; i++)
        {
            path_cm[i].reset();
            path_top[i].reset();
            w[i].reset();
        }
        for(int i = 0; i < INITAL_FRAME_NUMBER + 5; i++)
            path_bottom[i].reset();
        
        index_path_bottom = 0;
        index_path_cm = 0;
        index_path_top = 0;
        index_w = 0;

        printf("INFO: class: Bottle > resetBottle() - reset end\n");
    }

    int getIndexPath_cm()
    {
        return index_path_cm;
    }

    int getIndexPath_top()
    {
        return index_path_top;
    }

    int getIndexPath_bottom()
    {
        return index_path_bottom;
    }

private:
    int index_path_cm;
    int index_path_top;
    int index_path_bottom;
    int index_w;
};

Vector getInternalDivision(Vector path_top, Vector path_bottom)
{
    Vector ans;
    ans.x = (path_bottom.x - path_top.x) * (BOTTLE_HIGHT - CM_HIGHT) / BOTTLE_HIGHT + path_top.x;
    ans.y = (path_bottom.y - path_top.y) * (BOTTLE_HIGHT - CM_HIGHT) / BOTTLE_HIGHT + path_top.y;
    ans.t = path_top.t;
    return ans;
}

Vector getVelocity(Vector* path, int length)
{
    TGraph* trace_cm_x = new TGraph();
    TGraph* trace_cm_y = new TGraph();
    
    TF1* path_x_t = new TF1("x_t", "[0]*x+[1]", MIN_X, MAX_X);
    TF1* path_y_t = new TF1("y_t", "-4.9*x*x+[0]*x+[1]", MIN_Y, MAX_Y);

    double vy;
    double vx;

    Vector velocity;

    for(int i = 0; i < length; i++)
    {
        trace_cm_x->SetPoint(trace_cm_x->GetN(), path[i].t, path[i].x);
        trace_cm_y->SetPoint(trace_cm_y->GetN(), path[i].t, path[i].y);
        //printf("path_cm = %lf: (%lf, %lf)\n", path[i].t, path[i].x, path[i].y);
    }

    trace_cm_x->Fit(path_x_t, "M");
    trace_cm_y->Fit(path_y_t, "M");

    vx = path_x_t->GetParameter(0);
    vy = path_y_t->GetParameter(0);

    velocity.x = vx;
    velocity.y = vy;

    return velocity;
}

double innerProduct(Vector v1, Vector v2)
{
    return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}

Bottle errorCorrection(Bottle bottle)
{
    Bottle corrected_bottle;
    Vector corrected_vector;

    printf("-- error correction process start: top --\n");

    for(int i = 0; i < bottle.getIndexPath_top() - 2; i++)
    {
        corrected_vector.t = (bottle.path_top[i].t + bottle.path_top[i+1].t) / 2.0;
        corrected_vector.x = (bottle.path_top[i].x + bottle.path_top[i+1].x) / 2.0;
        corrected_vector.y = (bottle.path_top[i].y + bottle.path_top[i+1].y) / 2.0;
        //corrected_vector.print();
        corrected_bottle.setPath_top(corrected_vector);
    }

    printf("-- error correction process start: bottom --\n");
    for(int i = 0; i < bottle.getIndexPath_bottom() - 2; i++)
    {
        corrected_vector.t = (bottle.path_bottom[i].t + bottle.path_bottom[i+1].t) / 2.0;
        corrected_vector.x = (bottle.path_bottom[i].x + bottle.path_bottom[i+1].x) / 2.0;
        corrected_vector.y = (bottle.path_bottom[i].y + bottle.path_bottom[i+1].y) / 2.0;
        //corrected_vector.print();
        corrected_bottle.setPath_bottom(corrected_vector);
    }

    return corrected_bottle;    
}

FILE* fp_data; // file pointer for bottle path data file
FILE* fp_record;
Bottle bottle_raw;
Bottle bottle;
Vector top;  // temporary memory for bottle_top
Vector bottom; // temporary memory for bottle_bottom
Vector r_cm_1; // temporary memory for getting angular velocity
Vector r_cm_2; // temporary memory for getting angular velocity
Vector path_cm_temp; // temporary memory for recording path of CM
Vector angular_velocity; // tempoary memory for get angular velocity
Vector inital;
Vector last;

int n_line = 0;
int n_data = 0; // number of data by a bottle flip.
int experiment_no = 0;
int step = 0;
int bottom_step = 0;
int is_success = FALSE;

double x = 0.0;
double y = 0.0;
double dt = 0.0;
double t = 0.0;
double min_y = 0.0;
char graph_name[250];
char file_name[250] = {"rectangular_fail"};
char path[250];
double angular_displacement = 0.0;
//TMultiGraph* mg = new TMultiGraph();
TGraph* graph_top[MAX_EXPERIENCE_NUM];
TGraph* graph_bottom[MAX_EXPERIENCE_NUM];
TGraph* graph_cm[MAX_EXPERIENCE_NUM];
TGraph* graph_w[MAX_EXPERIENCE_NUM];
TMultiGraph* mg[MAX_EXPERIENCE_NUM];

void BottleFlip()
{
    TH1D* hist_angular_displacement = new TH1D("angular displacement", "angular displacement", 300, 3, 6);
    TH1D* hist_angular_velocity = new TH1D("angular velocity", "angular velocity", 50, 0, 50);
    TH2D* hist_inital_velocity = new TH2D("inital velocity", "inital velocity", 100, -2.0, 2.0, 100, -2.0, 2.0);

    for(int i = 0; i < MAX_EXPERIENCE_NUM; i++)
    {
        graph_top[i] = new TGraph();
        graph_bottom[i] = new TGraph();
        graph_cm[i] = new TGraph();
        graph_w[i] = new TGraph();
        mg[i] = new TMultiGraph();
    }
    // TGraph* graph_top = new TGraph();  // graph object for ploting bottle top
    // TGraph* graph_bottom = new TGraph();  // graph object for ploting bottle bottom
    // TGraph* graph_cm[experiment_no] = new TGraph();  // graph object for ploting center of mass

    auto c_w = new TCanvas("angular_velocity","angular_velocity",800, 600);
    auto c_path = new TCanvas("path_of_bottle","path_of_bottle",800, 600);
    auto c_angular_displacement = new TCanvas("angular_displacement","angular_displacement",800, 600);
    auto c_hist_w = new TCanvas("angular_velocity_distribution","angular_velocity_distribution",800, 600);
    auto c_inital_velocity = new TCanvas("inital_velocity","inital_velocity",800, 600);

    sprintf(path, "%s.csv", file_name);
    fp_data = fopen64(path, "r");
    sprintf(path, "%s_record.csv", file_name);
    fp_record = fopen64(path, "w");

    if(fp_data != nullptr && fp_record != nullptr)
    {
        while(!feof(fp_data)) 
        {
            n_line++;
            top.reset();
            bottom.reset();
            fscanf(fp_data, "%lf,%lf,%lf,%lf,%lf,%lf", &top.t, &top.x, &top.y, &bottom.t, &bottom.x, &bottom.y);

            // _inverse file only
            //top.x = (-1.0) * top.x;
            //bottom.x = (-1.0) * bottom.x;
            
            //printf("input: %lf,%lf,%lf,%lf,%lf,%lf\n", top.t, top.x, top.y, bottom.t, bottom.x, bottom.y);
            
            if(abs(top.t + 100000.0) < 0.1 || abs(top.x + 100000.0) < 0.1 || abs(top.y + 100000.0) < 0.1)
            {
                // TMultiGraph* mg = new TMultiGraph();
                // graph_top[experiment_no] = new TGraph();  // graph object for ploting bottle top
                // graph_bottom[experiment_no] = new TGraph();  // graph object for ploting bottle bottom
                // graph_cm[experiment_no] = new TGraph();  // graph object for ploting center of mass
                
                bottle = errorCorrection(bottle_raw);

                for(int i = 0; i < bottle.getIndexPath_top(); i++)
                    graph_top[experiment_no]->SetPoint(graph_top[experiment_no]->GetN(), bottle.path_top[i].x, bottle.path_top[i].y);
                
                for(int i = 0; i < bottle.getIndexPath_bottom(); i++)
                    graph_bottom[experiment_no]->SetPoint(graph_bottom[experiment_no]->GetN(), bottle.path_bottom[i].x, bottle.path_bottom[i].y);

                //printf("getIndexPath_bottom: %d\n", bottle.getIndexPath_bottom());
                for(int i = 0; i < bottle.getIndexPath_bottom(); i++)  // get path of cm from first (INITAL_FRAME_NUMBER) frame
                    bottle.setPath_cm(getInternalDivision(bottle.path_top[i], bottle.path_bottom[i]));
                   //bottle.path_cm[i] =;
                              
                bottle.r_cm_i = bottle.path_cm[bottle.getIndexPath_cm() - 1]; // set inital position of center of mass
                bottle.v_cm_i = getVelocity(bottle.path_cm, bottle.getIndexPath_bottom()); // get inital velocity of center of mass
                printf("INFO: v_cm_i = %lfi + %lfj\n", bottle.v_cm_i.x, bottle.v_cm_i.y);
                hist_inital_velocity -> Fill(bottle.v_cm_i.x, bottle.v_cm_i.y, bottle.v_cm_i.getSize());

                // TODO : set dt
                //dt = (bottle.path_top[bottle.getIndexPath_top()-1].t - bottle.path_top[0].t) / double(n_data);
                min_y = bottle.getMinY();
            
                for(int i = 0; i < bottle.getIndexPath_cm(); i++)
                    graph_cm[experiment_no] -> SetPoint(graph_cm[experiment_no]->GetN(), bottle.path_cm[i].x, bottle.path_cm[i].y);

                //printf("dt: %lf, n_data: %d\n", dt, n_data);
                step = bottle.getIndexPath_cm(); //reset step number
                do
                {
                    //dt = bottle.path_top[step+1].t - bottle.path_top[step].t;
                    //t = dt*step + bottle.r_cm_i.t;
                    t = bottle.path_top[step].t;
                    x = bottle.r_cm_i.x + bottle.v_cm_i.x * t;
                    y = bottle.r_cm_i.y + bottle.v_cm_i.y * t - g*t*t/2.0;

                    //printf("t: %lf, x: %lf, y: %lf\n", t, x, y);

                    path_cm_temp.x = x;
                    path_cm_temp.y = y;
                    path_cm_temp.t = t;
                    bottle.setPath_cm(path_cm_temp);

                    fprintf(fp_record, "%lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf, %lf\n", bottle.path_top[step].t, bottle.path_top[step].x, bottle.path_top[step].y, bottle.path_bottom[step].t, bottle.path_bottom[step].x, bottle.path_bottom[step].y, t, x, y);

                    graph_cm[experiment_no] -> SetPoint(graph_cm[experiment_no]->GetN(), x, y);
                    //printf("graph_cm[experiment_no]->GetN(): %d\n", int(graph_cm[experiment_no]->GetN()));
                    step++;

                }while(step < bottle.getIndexPath_top());
                fprintf(fp_record, "-100000, -100000, -100000, -100000, -100000, -100000, -100000, -100000, -100000\n");

                // Get angular velocity
                graph_w[experiment_no] -> SetPoint(graph_w[experiment_no] -> GetN(), (bottle.path_top[1].t - bottle.path_top[0].t)/2.0, 0.0);
                printf("INFO: step = %d\n", step);
                double t = 0.0;
                double w = 0.0;

                for(int i = 0; i < step - 1; i++)
                {
                    if(i > bottle.getIndexPath_top())
                    {
                        printf("over the data\n");
                        break;
                    }
                        
                    r_cm_1.x = bottle.path_top[i].x - bottle.path_cm[i].x;
                    r_cm_1.y = bottle.path_top[i].y - bottle.path_cm[i].y;
                    r_cm_1.t = bottle.path_top[i].t;

                    r_cm_2.x = bottle.path_top[i+1].x - bottle.path_cm[i+1].x;
                    r_cm_2.y = bottle.path_top[i+1].y - bottle.path_cm[i+1].y;
                    r_cm_2.t = bottle.path_top[i+1].t;

                    if(i == 0)
                    {
                        inital = r_cm_1;
                        printf("set inital direction vector\n");
                    }
                    else if(i == step - 2)
                    {
                        last = r_cm_2;
                        printf("set last direction vector\n");
                    }
                        

                    if(fabs(r_cm_1.t - r_cm_2.t) < MIN_TIME_INTERVER)
                    {
                        printf("ERROR!: Get angular velocity >> time interver < MIN_TIME_INTERVER\n");
                        continue;
                    }
                        

                    angular_velocity.z = acos(innerProduct(r_cm_1, r_cm_2)/(r_cm_1.getSize() * r_cm_2.getSize())) / (r_cm_2.t - r_cm_1.t);
                    angular_velocity.t = (r_cm_1.t + r_cm_2.t)/2.0;
                    bottle.setW(angular_velocity);
                    
                    graph_w[experiment_no] -> SetPoint(graph_w[experiment_no]->GetN(), (r_cm_1.t + r_cm_2.t)/2.0, angular_velocity.z);

                    //printf("angular velocity: %lf\n", angular_velocity.z);
                }

                //max_w = graph_w[experiment_no]->GetMaximum();
                graph_w[experiment_no]-> SetPoint(graph_w[experiment_no]->GetN(), (r_cm_1.t + r_cm_2.t)/2.0, 0);

                // calculus angular displacement
                angular_displacement = 2*M_PI - acos(innerProduct(inital, last) / (inital.getSize() * last.getSize()));
                hist_angular_displacement -> Fill(angular_displacement); 
                angular_displacement = angular_displacement / M_PI;

                // record inital angular velocity
                hist_angular_velocity -> Fill(bottle.w[0].z);

                // draw graph
                c_path -> cd();
                graph_top[experiment_no] -> SetTitle("path_top");
                graph_top[experiment_no] -> SetMarkerColor(2);
                graph_top[experiment_no] -> SetMarkerSize(1.0);
                graph_top[experiment_no] -> SetMarkerStyle(2);
                graph_top[experiment_no] -> SetLineColor(2);
               
                graph_bottom[experiment_no] -> SetTitle("path_bottom");
                graph_bottom[experiment_no] -> SetMarkerColor(3);
                graph_bottom[experiment_no] -> SetMarkerSize(1.0);
                graph_bottom[experiment_no] -> SetMarkerStyle(2);
                graph_bottom[experiment_no] -> SetLineColor(3);

                graph_cm[experiment_no] -> SetTitle("path_cm");
                graph_cm[experiment_no] -> SetMarkerColor(4);
                graph_cm[experiment_no] -> SetMarkerSize(1.0);
                graph_cm[experiment_no] -> SetMarkerStyle(2);
                graph_cm[experiment_no] -> SetLineColor(4);
               
                mg[experiment_no] -> Add(graph_top[experiment_no]);
                mg[experiment_no] -> Add(graph_bottom[experiment_no]);
                mg[experiment_no] -> Add(graph_cm[experiment_no]);
                mg[experiment_no] -> GetXaxis()->SetTitle("position: x[m]");
                mg[experiment_no] -> GetYaxis()->SetTitle("position: y[m]");
                mg[experiment_no] -> GetYaxis() -> SetRangeUser(MIN_Y, MAX_Y);
                mg[experiment_no] -> GetXaxis() -> SetLimits(MIN_X, MAX_X);

                mg[experiment_no] -> Draw("APL APL APL");
                c_path -> BuildLegend();
                
                if(is_success == TRUE)
                    sprintf(graph_name, "/home/lsh/Dropbox/GYPT/graph/success/%s_%d.png",file_name ,experiment_no);
                else
                    sprintf(graph_name, "/home/lsh/Dropbox/GYPT/graph/fail/%s_%d.png",file_name, experiment_no);

                
                c_path -> SaveAs(graph_name);

                c_w -> cd();
                sprintf(graph_name, "angular velocity - total angular displacement = %lf pi(rad)", angular_displacement);
                graph_w[experiment_no] -> SetTitle(graph_name);
                graph_w[experiment_no] -> SetMarkerColor(6);
                graph_w[experiment_no] -> SetMarkerSize(1.0);
                graph_w[experiment_no] -> SetMarkerStyle(2);
                graph_w[experiment_no] -> SetLineColor(6);
                graph_w[experiment_no] -> GetXaxis()->SetTitle("time [s]");
                graph_w[experiment_no] -> GetYaxis()->SetTitle("angular velocity [(rad)/s]");

                graph_w[experiment_no]->Draw("APL");

                if(is_success == TRUE)
                    sprintf(graph_name, "/home/lsh/Dropbox/GYPT/graph/angular_velocity/success/%s_angular_velocity_%d.png",file_name ,experiment_no);
                else
                    sprintf(graph_name, "/home/lsh/Dropbox/GYPT/graph/angular_velocity/fail/%s_angular_velocity_%d.png",file_name, experiment_no);
            
                c_w->SaveAs(graph_name);

                //reset
                n_data = 0;
                bottle.resetBottle();
                bottle_raw.resetBottle();
                c_path->Clear();
                c_w->Clear();
                bottom_step = 0;
                experiment_no++;

                printf("\n==================== %d - bottle end ====================\n", experiment_no);
                if(experiment_no -1 >= MAX_EXPERIENCE_NUM)
                {
                    printf("Error!!: number of experement is to large || experiment_no = %d\n", experiment_no);
                    exit(-1);
                }
                
            }
            else
            {
                top.t = TIME_INTERVER_ratio * top.t;
                bottom.t = TIME_INTERVER_ratio * bottom.t;
                bottle_raw.setPath_top(top);
                if(!(abs(bottom.x) < 0.0001) && !(abs(bottom.y) < 0.0001) && (bottom_step <= INITAL_FRAME_NUMBER))
                    bottle_raw.setPath_bottom(bottom);
                    bottom_step++;
                
                n_data++;
            }
        }

        c_angular_displacement -> cd();
        hist_angular_displacement -> SetTitle("angular displacement distribution");
        hist_angular_displacement -> GetXaxis() -> SetTitle("angular displacement (rad)");
        hist_angular_displacement -> GetYaxis() -> SetTitle("number");
        hist_angular_displacement -> Draw();

        sprintf(graph_name, "/home/lsh/Dropbox/GYPT/graph/angular_displacement_distribution_%s.png",file_name);   
        //c_angular_displacement -> SaveAs("/home/lsh/Dropbox/GYPT/graph/angular_displacement_distribution_round_success.png");
        c_angular_displacement -> SaveAs(graph_name);

        c_hist_w -> cd();
        hist_angular_velocity -> SetTitle("angular velocity distribution");
        hist_angular_velocity -> GetXaxis() -> SetTitle("angular velocity [rad/s]");
        hist_angular_velocity -> GetYaxis() -> SetTitle("number");
        hist_angular_velocity -> Draw();
        //c_hist_w -> SaveAs("/home/lsh/Dropbox/GYPT/graph/angular_velocity_distribution_round_success.png");
        sprintf(graph_name, "/home/lsh/Dropbox/GYPT/graph/angular_velocity_distribution_%s.png",file_name);
        c_hist_w -> SaveAs(graph_name);

        c_inital_velocity -> cd();
        hist_inital_velocity -> SetTitle("inital velocity(CM) distribution");
        hist_inital_velocity -> GetXaxis() -> SetTitle("inital x axis velocity [m/s]");
        hist_inital_velocity -> GetYaxis() -> SetTitle("inital y axis velocity [m/s]");
        hist_inital_velocity -> SetMarkerStyle(20);
        if(is_success == TRUE)
            hist_inital_velocity -> SetMarkerColor(kGreen);
        else
             hist_inital_velocity -> SetMarkerColor(kRed);
        hist_inital_velocity -> Draw();
        //c_inital_velocity -> SaveAs("/home/lsh/Dropbox/GYPT/graph/inital_velocity_distribution_round_success.png");
        sprintf(graph_name, "/home/lsh/Dropbox/GYPT/graph/inital_velocity_distribution_%s.png",file_name);
        c_inital_velocity -> SaveAs(graph_name);

        fclose(fp_data);
        fclose(fp_record);
    }
}