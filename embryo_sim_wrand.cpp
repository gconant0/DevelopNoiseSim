#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <math.h>
#include <algorithm>
#include "plot.h"
#include <map>


#define MAX_D 25
#define MAX_COLOR 65535

int num_cell_trans, trans_id;

using namespace std;

class Pow2 {
public:
    Pow2();
    Pow2(int d);
    int get_pow(int i);
    ~Pow2();
protected:
    int depth, *pow_array;
};

class Cell {
public:
    Cell(Pow2 *m);
    Cell(Cell* parent, int offs, Pow2 *m);
    int get_my_id()  {return(cell_id);};
    int get_my_global_id() {return(global_id);};
    int get_my_generation() {return(my_generation);};
    int get_my_offspring()  {return(my_offspring);};
    int get_my_lineage() {return(my_lineage);};
    int get_my_lineage_level() {return(my_lineage_level);};
    Cell* get_my_parent() {return(my_parent);};
    void set_child(Cell* new_child, int c_num);
    Cell* get_child(int child_num);
    void mark_failed () {failed_trans=true;};
    bool is_failed ()  {return(failed_trans);};
    void set_lineage(int lin_val, int lin_level_val);
    void set_coords(double y_level, double xspace, double rad, double border);
    
    void draw_cell(int **lineage_colors, int max_lineage,  int max_gen,  double rad);
    void draw_lines();
    double get_x()  {return(xcent);};
    double get_y()  {return(ycent);};
    
    
protected:
    int cell_id, global_id, my_lineage, my_lineage_level, my_generation, my_offspring;
    double xcent, ycent;
    bool failed_trans;
    Cell *my_parent, *child[2];
    Pow2 *my_pow2;
    
    void recurse_coords(double y_level, double xspace, double rad, double border);
    void recurse_draw(int **lineage_colors, int max_lineage, int max_gen,  double rad);
    void recurse_lines();
};


void divide (Cell *my_cell, int max_d, Pow2 *my_pow2);
void delete_lineage (Cell *my_cell);
void fill_final_gen(Cell *my_cell, Cell **final_gen);
void assign_lineages(Cell *my_cell, int *trans_points, int *lineage_levels,int num_trans, int gens_per_tran, Pow2 *my_pows);
void assign_lineages_rand(Cell *my_cell, int *trans_points, int *lineage_levels, int num_trans, int gens_per_tran, Pow2 *my_pows, int *rand_vals, int cut, bool robust_trans);
int make_diagram (Cell *root_cell, string outfile, Pow2 *my_pow2, int num_trans, int total_depth);

int main(int argc, char** argv)
{
    int i, j, num_trans, gens_per_tran, num_cells, total_cells, num_generations, l, num_sims, num_locs=0, spp_pos, *factor_locs, lpos,
    *trans_points, num_lineages, *lineage_cnt, int_error, *prot_vals, val, *sorted_vals, loc_cut, min_val, *lineage_vals;
    double vect_len, *base_vect, error, *new_vect, new_len, dist;
    string runcmd, model, output_dir, outputfile, line, newline, cutnum, spp_name;
    bool use_given_cut=false, robust_trans=false;
    std::map<std::string, int> factloc_map;
    ifstream fin;
    Pow2 *my_pow2;
    Cell *base_cell, *new_cell, *parent, **final_gen, *base_rand;
    
    if (argc> 3) {
        std::stringstream ss(argv[1]);
        ss >> num_trans;
        
        std::stringstream ss2(argv[2]);
        ss2 >> gens_per_tran;
        
        model=argv[3];
        
        output_dir=model;
        output_dir = output_dir.substr(0, output_dir.length()-4);
        
        output_dir=output_dir + "_output/trajectories/";
        
        cout<<"Output directory will be "<<output_dir<<endl;
        
        if ((argv[4][1] == 'n') || (argv[4][1] == 'N')) {
            cutnum=argv[4];
            cutnum = cutnum.substr(3, cutnum.length()-3);
            std::stringstream ss5(cutnum);
            ss5 >> min_val;
            use_given_cut=true;
            cout<<"Factor cutoff value (from "<<cutnum<<") is "<<min_val<<endl;
        }
        else {
            std::stringstream ss3(argv[4]);
            ss3 >> error;
            error=fabs(error);
            int_error=(int)error;
            error = error - (double)int_error;
            cout<<"Error cutoff is "<<error<<endl;
        }
        
        
        if (argc>5) robust_trans=true;
        
        fin.open(model.c_str());
        if (!(fin.fail())) {
            spp_pos=0;
            while(!fin.eof()) {
                getline(fin,newline);
                if (!(fin.eof())) {
                    size_t pos = newline.find("<Species>");

                    if (pos != std::string::npos) {
                        getline(fin,newline);
                        
                        size_t pos2 = newline.find("<Id>");

                        if (pos2 != std::string::npos) {
                            newline.erase(pos2, 4);
                        }
                        
                        size_t pos3 = newline.find("</Id>");

                        if (pos3 != std::string::npos) {
                            newline.erase(pos3, 5);
                        }
                        
                        std::stringstream ss(newline);
                        ss >> spp_name;
                        
                        //cout<<"Species # "<<spp_pos<<" is "<<spp_name<<endl;
                        
                        if (spp_name == "TF") {factloc_map[spp_name]=spp_pos; num_locs++;}
                        if (spp_name == "TF1") {factloc_map[spp_name]=spp_pos; num_locs++;}
                        if (spp_name == "TF2") {factloc_map[spp_name]=spp_pos; num_locs++;}
                        
                        spp_pos++;
                    }
                }
            }
            fin.close();
        }
        else {
            cerr<<"ERROR: Cannot open model file "<<model<<". Exiting"<<endl;
            return(-1);
        }
        
        cout<<"There are "<<num_locs<<" locations of the transition factor"<<endl;
        
        factor_locs=new int [num_locs];
        
        i=0;
        for (auto it = factloc_map.begin(); it != factloc_map.end(); ++it) {
            factor_locs[i]=it->second;
            factor_locs[i]++;
            cout<<"Placing factor # "<<i<<" at "<<factor_locs[i]<<endl;
            i++;
        }
        
        num_trans=abs(num_trans);
        gens_per_tran=abs(gens_per_tran);
        
        trans_points=new int[num_trans];
        
        trans_points[0]=gens_per_tran;
        for(i=1; i<num_trans; i++)
            trans_points[i]=trans_points[i-1]+gens_per_tran;
        
        for(i=0; i<num_trans; i++) cout<<i<<": "<<trans_points[i]<<endl;
        
        my_pow2=new Pow2(MAX_D);
        
        base_cell= new Cell(my_pow2);
        
        parent=base_cell;
        num_generations=(num_trans+1)*gens_per_tran;
        
        num_lineages=my_pow2->get_pow(num_trans);
        lineage_cnt=new int[num_lineages];
        base_vect=new double[num_lineages];
        
        for(i=0; i<num_lineages; i++) lineage_cnt[i]=0;
        
        num_cells=my_pow2->get_pow(num_generations-1);
        cout<<"Expecting "<<num_cells<<" cells in the final generation ("<<num_generations<<")"<<" There are "<<num_lineages<<" lineages"<<endl;
        
        total_cells=2*num_cells-1;
        cout<<"System has "<<total_cells<<" total cells"<<endl;
        
        lineage_vals=new int[num_generations];
        
        l=0;
        for(i=0; i<num_generations; i++) {
            if (i == trans_points[l])
                l++;
            lineage_vals[i]=l;
        }
        
        
        num_cell_trans=0;
        divide(base_cell, num_generations, my_pow2);
            // for(i=0; i<MAX_D; i++) cout<<"i: "<<i<<" POW2: "<<my_pow2->get_pow(i)<<endl;
        assign_lineages(base_cell->get_child(0), trans_points, lineage_vals, num_trans, gens_per_tran, my_pow2);
        assign_lineages(base_cell->get_child(1), trans_points, lineage_vals, num_trans, gens_per_tran, my_pow2);
        
        std::stringstream ss4;
        num_sims=num_cell_trans;
        
        if (robust_trans==true) num_sims=total_cells-1;
        
        ss4<<"tau_leaping_exp_adapt_serial -m "<<model<<" -t 3000 -r "<<num_sims<<" -i 100 --force --keep-trajectories";
        runcmd=ss4.str();
        cout<<"Running "<<runcmd<<endl;
        
        int ret = std::system(runcmd.c_str());

        if (ret != 0) {
            std::cerr << "Command failed with return code " << ret << std::endl;
            return(-1);
        }
        
        prot_vals=new int[num_cell_trans];
        sorted_vals=new int[num_cell_trans];
        
        for(i=0; i<num_cell_trans; i++) {
            prot_vals[i]=0;
            std::stringstream ss5;
            ss5<<output_dir<<"trajectory"<<i<<".txt";
            outputfile=ss5.str();
            fin.open(outputfile.c_str());
            if (!(fin.fail())) {
                while(!fin.eof()) {
                    getline(fin,newline);
                    if (!(fin.eof())) line=newline;
                }
                fin.close();
                //cout<<"Last data: "<<line<<endl;
                std::stringstream ss6(line);
                lpos=0;
                j=0;
                while((j<=factor_locs[lpos]) && (lpos<num_locs)) {
                    ss6>>val;

                    //for(j=0; j<9; j++) ss6>>val;
                    if (j == factor_locs[lpos]) {
                        //cout<<"At j = "<<j<<" for lpos = "<<lpos<<" meaning "<<factor_locs[lpos]<<" Read: "<<val<<endl;
                        prot_vals[i]+=val;
                        lpos++;
                    }
                    j++;
                }
                //cout<<"i "<<i<<" v="<<val<<endl;
            }
            
        }
        
        if (use_given_cut == false) {
            for(i=0; i<num_cell_trans; i++) sorted_vals[i]=prot_vals[i];
            std::sort(sorted_vals, sorted_vals + num_cell_trans);
            
            loc_cut=(int)((1.0-error)*(double)num_cell_trans);
            
            min_val=sorted_vals[loc_cut];
            cout<<"TF val cutoff of error "<<error<<" is at "<<loc_cut<<" which is "<<min_val<<"( "<<sorted_vals[loc_cut-1]<<", "<<sorted_vals[loc_cut+1]<<")"<<endl;
            
        }
      
        
        
        final_gen=new Cell*[num_cells];
        fill_final_gen(base_cell, final_gen);
        
        make_diagram (base_cell, "perfectcell.ps", my_pow2, num_trans, num_generations);
        
        
        //cout<<"CellID\tParID\tLineageID\tLineageLevel"<<endl;
        for(i=0; i<num_cells; i++) {
            //cout<<final_gen[i]->get_my_id()<<"\t"<<final_gen[i]->get_my_parent()->get_my_id()<<"\t"<<final_gen[i]->get_my_lineage()<<"\t"<<final_gen[i]->get_my_lineage_level()<<endl;
            if (final_gen[i]->get_my_lineage_level()==num_trans)
                lineage_cnt[final_gen[i]->get_my_lineage()]++;
        }
        
        vect_len=0;
        cout<<"Lineage\tCnt"<<endl;
        for(i=0; i<num_lineages; i++) {
            //cout<<i<<"\t"<<lineage_cnt[i]<<endl;
            vect_len+=((double)lineage_cnt[i]*(double)lineage_cnt[i]);
        }
        
        
        vect_len=sqrt(vect_len);
        cout<<"Length of final system: "<<vect_len<<endl;
       
        for(i=0; i<num_lineages; i++) base_vect[i]=(double)lineage_cnt[i]/vect_len;
        
        base_rand=new Cell(my_pow2);
        divide(base_rand, num_generations, my_pow2);
        trans_id=0;
        assign_lineages_rand(base_rand->get_child(0), trans_points, lineage_vals, num_trans, gens_per_tran, my_pow2, prot_vals, min_val, robust_trans);
        assign_lineages_rand(base_rand->get_child(1), trans_points, lineage_vals, num_trans, gens_per_tran, my_pow2, prot_vals, min_val, robust_trans);
        
        fill_final_gen(base_rand, final_gen);
        
        make_diagram (base_rand, "randcell.ps", my_pow2, num_trans, num_generations);
        
        
        //for(i=0; i<num_cells; i++)
           // cout<<final_gen[i]->get_my_id()<<"\t"<<final_gen[i]->get_my_parent()->get_my_id()<<"\t"<<final_gen[i]->get_my_lineage()<<"\t"<<final_gen[i]->get_my_lineage_level()<<endl;
        
        new_vect=new double[num_lineages];
        
        for(i=0; i<num_lineages; i++) {
            lineage_cnt[i]=0;
            new_vect[i]=0.0;
        }
        
        for(i=0; i<num_cells; i++) {
            if (final_gen[i]->get_my_lineage_level()==num_trans)
                lineage_cnt[final_gen[i]->get_my_lineage()]++;
        }
        
        new_len=0;
        for(i=0; i<num_lineages; i++)
            new_len+=((double)lineage_cnt[i]*(double)lineage_cnt[i]);
        new_len=sqrt(new_len);
        cout<<"Length of rand system: "<<new_len<<endl;
        for(i=0; i<num_lineages; i++) new_vect[i]=(double)lineage_cnt[i]/new_len;
       
        dist=0;
        
        for(i=0; i<num_lineages; i++) dist+= ((base_vect[i]-new_vect[i])*(base_vect[i]-new_vect[i]));
        
        dist=sqrt(dist);
        
        cout<<"Dist between orig and rand: "<<dist<<endl;
        
        delete[] lineage_cnt;
        delete[] base_vect;
        delete[] final_gen;
        delete_lineage(base_cell);
        delete_lineage(base_rand);
        delete my_pow2;
    }
    else {
        cerr<<"Usage: embryo_sim #CellTransitions #GenerationsPerTrans <ModelFile> Error"<<endl;
    }


}


void divide (Cell *my_cell, int max_d, Pow2 *my_pow2)
{
    Cell *child1, *child2;
    
    child1=new Cell(my_cell, 0, my_pow2);
    child2=new Cell(my_cell, 1, my_pow2);
    
    my_cell->set_child(child1, 0);
    my_cell->set_child(child2, 1);
    
    if (child1->get_my_generation() < (max_d-1)) divide(child1, max_d, my_pow2);
    if (child2->get_my_generation() < (max_d-1)) divide(child2, max_d, my_pow2);
}


void fill_final_gen(Cell *my_cell, Cell **final_gen)
{
    if (my_cell->get_child(0) !=0) fill_final_gen(my_cell->get_child(0), final_gen);
    else final_gen[my_cell->get_my_id()]=my_cell;
    if (my_cell->get_child(1) !=0) fill_final_gen(my_cell->get_child(1), final_gen);
    else final_gen[my_cell->get_my_id()]=my_cell;
}


void assign_lineages(Cell *my_cell, int *trans_points,int *lineage_levels, int num_trans, int gens_per_tran, Pow2 *my_pows)
{
    int genloc=0, new_level, new_lineage_level, bit_mask, my_bit;
    
    if (my_cell->get_my_generation()> trans_points[0]) {
        while((trans_points[genloc] < my_cell->get_my_generation()) && (genloc<num_trans)) genloc++;
    }
  
    
    if (lineage_levels[my_cell->get_my_generation()] == (lineage_levels[my_cell->get_my_parent()->get_my_generation()]+1)) {
   //if (my_cell->get_my_generation() == trans_points[genloc]) {
        
        
        new_lineage_level=my_cell->get_my_parent()->get_my_lineage_level()+1;
        new_level=my_cell->get_my_parent()->get_my_lineage()*2;
        
        if (gens_per_tran ==1) {
            if (my_cell->get_my_parent()->get_child(1)==my_cell)
                new_level++;
        }
        else {
            bit_mask=(my_pows->get_pow(my_cell->get_my_parent()->get_my_lineage_level()) & my_cell->get_my_parent()->get_my_id());
            my_bit = bit_mask >> (my_cell->get_my_parent()->get_my_lineage_level());
            new_level+=my_bit;
        }
        //if (my_cell == my_cell->get_my_parent()->get_child(1)) new_level++;
        
        my_cell->set_lineage(new_level, new_lineage_level);
        num_cell_trans++;
        
    }
    else
        my_cell->set_lineage(my_cell->get_my_parent()->get_my_lineage(), my_cell->get_my_parent()->get_my_lineage_level());
    
    if (my_cell->get_child(0) != 0) assign_lineages(my_cell->get_child(0), trans_points, lineage_levels, num_trans, gens_per_tran, my_pows);
    if (my_cell->get_child(1) != 0) assign_lineages(my_cell->get_child(1), trans_points, lineage_levels, num_trans, gens_per_tran, my_pows);
    
}


void assign_lineages_rand(Cell *my_cell, int *trans_points, int *lineage_levels, int num_trans, int gens_per_tran, Pow2 *my_pows, int *rand_vals, int cut, bool robust_trans)
{
    int genloc=0, new_level, new_lineage_level, bit_mask, my_bit;
    bool should_trans=false, failed_trans=false;
    
    
    my_cell->set_lineage(my_cell->get_my_parent()->get_my_lineage(), my_cell->get_my_parent()->get_my_lineage_level());
    
  
    if (robust_trans==false) {
        if (lineage_levels[my_cell->get_my_generation()] == (lineage_levels[my_cell->get_my_parent()->get_my_generation()]+1)) should_trans=true;
    }
    else {
        if (lineage_levels[my_cell->get_my_generation()] > (my_cell->get_my_parent()->get_my_lineage_level())) {
            should_trans=true;
        }
    }
    
    if (should_trans) {
        if (rand_vals[trans_id]<cut) failed_trans=true;
    }

    if ((should_trans==true) && (failed_trans==false)) {
        new_lineage_level=my_cell->get_my_parent()->get_my_lineage_level()+1;
        new_level=my_cell->get_my_parent()->get_my_lineage()*2;
        
        if ((gens_per_tran ==1) || (robust_trans == true)) {
            if (my_cell->get_my_parent()->get_child(1)==my_cell)
                new_level++;
        }
        else {
            
            bit_mask=(my_pows->get_pow(my_cell->get_my_parent()->get_my_lineage_level()) & my_cell->get_my_parent()->get_my_id());
            my_bit = bit_mask >> (my_cell->get_my_parent()->get_my_lineage_level());
            new_level+=my_bit;
        }
        
        //if (my_cell == my_cell->get_my_parent()->get_child(1)) new_level++;
        
        my_cell->set_lineage(new_level, new_lineage_level);
        
    }
    
    
    if (failed_trans==true) my_cell->mark_failed ();
    if (should_trans==true) trans_id++;
    
    
    if (my_cell->get_child(0) != 0) assign_lineages_rand(my_cell->get_child(0), trans_points, lineage_levels,  num_trans, gens_per_tran, my_pows, rand_vals, cut, robust_trans);
    if (my_cell->get_child(1) != 0) assign_lineages_rand(my_cell->get_child(1), trans_points, lineage_levels, num_trans, gens_per_tran, my_pows, rand_vals, cut, robust_trans);
}



void delete_lineage (Cell *my_cell)
{
    if (my_cell->get_child(0) !=0) delete_lineage(my_cell->get_child(0));
    if (my_cell->get_child(1) !=0) delete_lineage(my_cell->get_child(1));
    delete my_cell;
    
}



Pow2::Pow2()
{
    int i;
    depth=31;
    pow_array=new int[depth];
    
    pow_array[0]=1;
    for(i=1; i<depth; i++)
        pow_array[i]=2*pow_array[i-1];
    
}



Pow2::Pow2(int d)
{
    int i;
    
    d=abs(d);
    
    depth=(d%31);
    pow_array=new int[depth];
    
    pow_array[0]=1;
    for(i=1; i<depth; i++)
        pow_array[i]=2*pow_array[i-1];
    
}



int Pow2::get_pow(int i)
{
    return(pow_array[(abs(i)%depth)]);
}



Pow2::~Pow2()
{
    delete[] pow_array;
}


Cell::Cell(Pow2 *m)
{
    my_pow2=m;
    cell_id=0;
    global_id=0;
    my_lineage=0;
    my_lineage_level=0;
    my_parent=0;
    my_generation=0;
    my_offspring=0;
    child[0]=0;
    child[1]=0;
    failed_trans=false;
    
}



Cell::Cell(Cell* parent, int offs, Pow2 *m)
{
    my_pow2=m;
    my_offspring=(offs%2);
    global_id=parent->get_my_global_id()+my_offspring;
    child[0]=0;
    child[1]=0;
    failed_trans=false;
    
    if (parent !=0)
    {
        my_generation=parent->get_my_generation()+1;
        cell_id=(parent->get_my_id()*2)+my_offspring;
        my_parent=parent;
    }
    else {
        cell_id=0;
        my_lineage=0;
        my_parent=0;
        my_generation=0;
    }
}

void Cell::set_lineage(int lin_val, int lin_level_val)
{
    my_lineage=lin_val;
    my_lineage_level=lin_level_val;
}


void Cell::set_child(Cell* new_child, int c_num)
{
    int child_num;
    child_num=(abs(c_num)%2);
    child[child_num]=new_child;
}


Cell* Cell::get_child(int child_num)
{
    return(child[(abs(child_num)%2)]);
}
#define MAX_CELL_ROW 64

void Cell::set_coords(double y_level, double spacex, double rad, double border)
{
    xcent=spacex/2.0;
    ycent=border;
    
    child[0]->recurse_coords(y_level, spacex, rad, border);
    child[1]->recurse_coords(y_level, spacex, rad, border);
    
}

void Cell::recurse_coords(double y_level, double spacex, double rad, double border)
{
    int my_off, layer_size;
    if (my_pow2->get_pow(my_generation)<=MAX_CELL_ROW)
        //ycent=my_parent->get_y()+y_level;
        ycent=my_generation*y_level+border;
    else {
        layer_size=my_pow2->get_pow(my_generation)/MAX_CELL_ROW;
        my_off=cell_id%layer_size;
        ycent=my_generation*y_level+border+ (double)my_off*2.1*rad;
        //ycent=my_parent->get_y()+y_level + (double)my_off+2.1*rad;
    }
    xcent= (spacex-2.0*border)* ((double)(cell_id)/(double)my_pow2->get_pow(my_generation)+0.5*(1.0/(double)my_pow2->get_pow(my_generation)));
    
    if (child[0]!=0)
        child[0]->recurse_coords(y_level, spacex, rad, border);
    if (child[1]!=0)
        child[1]->recurse_coords(y_level, spacex, rad, border);
}

void Cell::draw_lines()
{
    recurse_lines();
}

void Cell::draw_cell(int **lineage_colors, int max_lineage, int max_gen, double rad)
{
    recurse_draw(lineage_colors, max_lineage, max_gen, rad);
}

void Cell::recurse_lines()
{
    if (child[0] !=0) {
        pl_fline(xcent, ycent, child[0]->get_x(), child[0]->get_y());
        child[0]->recurse_lines();
    }
    
    if (child[1] !=0) {
        pl_fline(xcent, ycent, child[1]->get_x(), child[1]->get_y());
        child[1]->recurse_lines();
    }
    
}

void Cell::recurse_draw(int **lineage_colors, int max_lineage, int max_gen, double rad)
{
    int my_loc, b_int, g_int, r_int, my_level, max_level;
    double my_g, my_b, my_r;
    double text_width;
    string label;
    stringstream ss;
   //cout<<"C: "<<cell_id<<" MG: "<<my_generation<<" x= "<<xcent<<" y= "<<ycent<<" MAX GNE: "<<max_gen<<" MAX LIN: "<<max_lineage<<"My Lin: "<<my_lineage<<" MY LL: "<<my_lineage_level<<endl;
    my_level=my_pow2->get_pow(my_lineage_level)-1;
    max_level =my_pow2->get_pow(max_lineage)-1;
    
    
    if (my_lineage_level !=0) {
        my_b = (((double)my_lineage)/(double)(my_pow2->get_pow(my_lineage_level)-1.0));
        my_g = 1.0 - my_b;
    }
    else {
        my_b=my_g=0.5;
    }
    
    b_int = int(my_b*lineage_colors[my_level][1]);
    g_int = int(my_g*lineage_colors[my_level][1]);
    my_r = 1.0-((double)(my_lineage_level)/(double)(max_lineage));
    r_int=int(my_r*0.6*lineage_colors[max_level][1]);
    
    //cout<<"C: "<<cell_id<<" MAXL: "<<max_lineage<<" MG: "<<my_generation<<" My LL: "<<my_lineage_level<<" MY L: "<<my_level<<"My B: "<<my_b<<" = "<<b_int<<" My g: "<<my_g<<" = "<<g_int<<" My R: "<<my_r<<" = "<<r_int<<endl;
    
    if (failed_trans == true)
        pl_fillcolor(MAX_COLOR, 0, 0);
    else {
        if (my_lineage_level==max_lineage) {
           
            
            
            pl_fillcolor(r_int, g_int, b_int);
        }
        else {
            if (my_generation != max_gen) {
                my_loc = (int)(((my_lineage+0.50)/(double)my_pow2->get_pow(my_lineage_level))*(double)my_pow2->get_pow(max_lineage));
                pl_fillcolor(r_int, g_int, b_int);
            }
            else pl_fillcolor(MAX_COLOR, MAX_COLOR, MAX_COLOR);
        }
    }
    pl_fcircle (xcent, ycent, rad);
    
    ss<<"ID: "<<cell_id<<" L: "<<my_lineage<<" LL: "<<my_lineage_level;
    label=ss.str();
    text_width = pl_flabelwidth (label.c_str());
    pl_fmove (xcent+1.2*rad+0.6*text_width,ycent-2.0*rad);
    pl_alabel('c', 'c', label.c_str());
    
    if (child[0]!=0)
        child[0]->recurse_draw(lineage_colors, max_lineage, max_gen, rad);
    if (child[1]!=0)
        child[1]->recurse_draw(lineage_colors, max_lineage, max_gen, rad);
    
}




int make_diagram (Cell *root_cell, string file, Pow2 *my_pow2, int num_trans, int total_depth)
{
    int *y_levels, i, num_lineages, **lineage_colors, thandle;
    double border=10, y_size, spacex, spacey, circle_radx, circle_rady, circle_rad;
    FILE *outfile;

   
    outfile=fopen(file.c_str(), "w");
    
    spacex=1000.0;
    spacey=1000.0;
    
    
    y_size=(0.4*spacey-2.0*border)/(double)total_depth;
    
    if (total_depth > 8) y_size = y_size *1.2;
   
    cout<<"For "<<total_depth<<" generations the layer size is "<<y_size<<endl;
    
    circle_rady=0.3*y_size;
   
    if (my_pow2->get_pow(total_depth) <= MAX_CELL_ROW)
        circle_radx= 0.46*(spacex/(double)my_pow2->get_pow(total_depth));
    else
        circle_radx= 0.46*(spacex/(double)MAX_CELL_ROW);
                           
    if (circle_radx<circle_rady) circle_rad=circle_radx;
    else circle_rad=circle_rady;
    
    root_cell->set_coords(y_size, spacex, circle_rad, border);
    
    num_lineages=my_pow2->get_pow(num_trans);
    
    lineage_colors=new int*[num_lineages];
    
    for(i=0; i<num_lineages; i++) {
        lineage_colors[i]=new int [3];
        lineage_colors[i][1]=((double)(i+1)/(double)(num_lineages))*MAX_COLOR;
        lineage_colors[i][0]=0;
        lineage_colors[i][2]=(1.0-((double)(i+1)/(double)(num_lineages)))*MAX_COLOR;
       // cout<<i<<": R="<<lineage_colors[i][0]<<" G="<<lineage_colors[i][1]<<" B="<<lineage_colors[i][1]<<endl;
    }
    
    for(i=0; i<num_trans; i++) {
        //lineage_colors[i][0]=(int)((double)MAX_COLOR*(double)(num_trans-i)/(double)num_trans);
        lineage_colors[i][0]=(int)((double)MAX_COLOR*((double)(i+1)/(double)num_trans));
    
    }
 
    /* create a Postscript Plotter that writes to standard output */
    if ((thandle = pl_newpl ("ps", stdin, outfile, stderr)) < 0)
    {
        fprintf (stderr, "Couldn't create Plotter\n");
        return 1;
    }
    pl_selectpl (thandle);       /* select the Plotter for use */
    
    if (pl_openpl () < 0)       /* open Plotter */
    {
        fprintf (stderr, "Couldn't open Plotter\n");
        return 1;
    }
    pl_fspace (0.0, 0.0, spacex, spacey); /* specify user coor system */
    pl_flinewidth (0.1);       /* line thickness in user coordinates */
    pl_pencolorname ("black");    /* path will be drawn in red */
    pl_erase ();
    
    
    pl_fontname("Helvetica-Oblique");
    pl_ftextangle (0);         /* text inclination angle (degrees) */
    pl_ffontsize (8);
    
    pl_filltype(0);
    root_cell->draw_lines();
    pl_filltype(1);
    root_cell->draw_cell(lineage_colors, num_trans, total_depth-1, circle_rad);
    
    if (pl_closepl () < 0)      /* close Plotter */
    {
        fprintf (stderr, "Couldn't close Plotter\n");
        return 1;
    }
    
    pl_selectpl (0);            /* select default Plotter */
    if (pl_deletepl (thandle) < 0) /* delete Plotter we used */
    {
        fprintf (stderr, "Couldn't delete Plotter\n");
        return 1;
    }

    
    return(0);
}
