#include "Integrator.h"

Integrator::Integrator(double dt, double m, int numpart, int dim, Atom *atoms, double b, int Nc, double sig, double eps, double T_bath):
    dt(dt),
    m(m),
    numpart(numpart),
    dim(dim),
    atoms(atoms),
    b(b),
    Nc(Nc),
    sig(sig),
    eps(eps),
    T_bath(T_bath)
{
}

void Integrator::integrate()
{
    int i,t,j,k,l,g,h,o,p,q,i2,j2,k2;
    int n = 300;
    vec v (dim);
    vec r (dim);
    vec r2 (dim);
    vec d (dim);
    vec F (dim);
    vec F_old (dim);
    double pot = 0.0;
    double E_k = 0.0;
    double E_tot;
    double T;
    vec v_half (dim);
    vec r_new (dim);
    vec v_new (dim);
    vec r_min (dim);
    double d_length;
    double L = b*(Nc);
    double V = L*L*L;
    double Fr_sum = 0.0;
    double P;
    vec r_init (dim);
    vec r_disp (dim);
    vec r_disp_min (dim);
    double r_disp_2;
    double r_sum = 0.0;
    double mean_disp;
    double D;
    int bins = 50;
    vec g_r (bins);
    int dist;
    double volume;
    double pi = 4*atan(1.0);
    double tau = 10*dt;
    double gamma;
    double tall;

    double r_cut = 3.0;
    int cells_x = L/r_cut;
    int nx,ny,nz;
    std::vector<Atom*> celle;
    std::vector<Atom*> celle2;

    Cell_container Con(cells_x);


    for(i=0;i<numpart;i++){
        r = atoms[i].getPosition();
        nx = int (r(0)/r_cut);
        ny = int (r(1)/r_cut);
        nz = int (r(2)/r_cut);

        if(nx>cells_x-1){
            nx = cells_x-1;
        }
        if(ny>cells_x-1){
            ny = cells_x-1;
        }
        if(nz>cells_x-1){
            nz = cells_x-1;
        }

        Con.container(nx,ny,nz,&atoms[i]);
    }



    ofstream myfile;
    myfile.open ("../Prosjekt 1/thermostat/andersen_gr_300.xyz");


    myfile << "Andersen thermostat.\n";

    for(i=0;i<bins;i++){
        g_r(i) = 0.0;
    }


    for(i=0;i<cells_x;i++){  // Finn kraft første gang
        for(j=0;j<cells_x;j++){
            for(k=0;k<cells_x;k++){
                celle = Con.getContainer(i,j,k);

                for(l=0;l<celle.size();l++){
                    r = celle[l]->getPosition();

                    for(g=0;g<dim;g++){
                        F(g) = 0;
                    }

                    for(g=0;g<celle.size();g++){  // egen celle
                        if(g!=l){
                            r2 = celle[g]->getPosition();
                            d = r - r2;
                            d_length = sqrt(d(0)*d(0) + d(1)*d(1) + d(2)*d(2));

                            pot = pot + 4*(pow (1/d_length,12) - pow (1/d_length,6));

                            for(h=0;h<dim;h++){
                                F(h) = F(h) + 24*(2*pow (1/d_length,12) - pow (1/d_length,6))*d(h)/(d_length*d_length);
                                Fr_sum = Fr_sum + F(h)*d(h);
                            }
                        }
                    }

                    for(g=-1;g<2;g++){   // naboceller
                        for(h=-1;h<2;h++){
                            for(o=-1;o<2;o++){
                                if(!(g==0 && h==0 && o==0)){

                                    i2 = i+g;
                                    j2 = j+h;
                                    k2 = k+o;

                                    if(i+g < 0){
                                        i2 = cells_x-1;
                                    }
                                    if(j+h < 0){
                                        j2 = cells_x-1;
                                    }
                                    if(k+o < 0){
                                        k2 = cells_x-1;
                                    }

                                    if(i+g == cells_x){
                                        i2 = 0;
                                    }
                                    if(j+h == cells_x){
                                        j2 = 0;
                                    }
                                    if(k+o == cells_x){
                                        k2 = 0;
                                    }

                                    celle2 = Con.getContainer(i2,j2,k2);

                                    for(p=0;p<celle2.size();p++){
                                        r2 = celle2[p]->getPosition();

                                        d = r-r2;

                                        for(q=0;q<dim;q++){   // minimum image convention
                                            if(fabs(d(q)-L) < fabs(d(q))){
                                               if(fabs(d(q)-L) < fabs(d(q)+L)){
                                                   r_min(q) = d(q)-L;
                                               }
                                            }
                                            if(fabs(d(q)) < fabs(d(q)+L)){
                                               if(fabs(d(q)) < fabs(d(q)-L)){
                                                  r_min(q) = d(q);
                                               }
                                            }
                                            if(fabs(d(q)+L) < fabs(d(q)-L)){
                                               if(fabs(d(q)+L) < fabs(d(q))){
                                                  r_min(q) = d(q)+L;
                                               }
                                            }
                                          }

                                        d_length = sqrt(r_min(0)*r_min(0) + r_min(1)*r_min(1) + r_min(2)*r_min(2));

                                        pot = pot + 4*(pow (1/d_length,12) - pow (1/d_length,6));

                                        for(q=0;q<dim;q++){
                                            F(q) = F(q) + 24*(2*pow (1/d_length,12) - pow (1/d_length,6))*r_min(q)/(d_length*d_length);
                                            Fr_sum = Fr_sum + F(q)*r_min(q);
                                        }
                                    }
                                }
                            }
                        }
                    }

                    celle[l]->setForce(F);
                }

            }
        }
    }




    for(t=0;t<n;t++){  // integrerer

        cout << t << endl;

//        stringstream fileNameStream;

//        fileNameStream << "../oppgi/oppgi004_" << t << ".xyz";
//        myfile.open (fileNameStream.str().c_str());
//        myfile << numpart << endl;
//        myfile << "Argon atoms in fcc lattice.\n";


        for(i=0; i<numpart;i++){  // Finner nye posisjoner
//            myfile << "Ar";

            r = atoms[i].getPosition();
            v = atoms[i].getVelocity();
            F = atoms[i].getForce();

            E_k = E_k + 1.0/2.0*(v(0)*v(0) + v(1)*v(1) + v(2)*v(2));


            r_init = atoms[i].getPosInit();

            for(j=0;j<dim;j++){
                r_disp(j) = r(j) - r_init(j);
            }

            for(q=0;q<dim;q++){   // minimum image convention
                if(fabs(r_disp(q)-L) < fabs(r_disp(q))){
                   if(fabs(r_disp(q)-L) < fabs(r_disp(q)+L)){
                       r_disp_min(q) = r_disp(q)-L;
                   }
                }
                if(fabs(r_disp(q)) < fabs(r_disp(q)+L)){
                   if(fabs(r_disp(q)) < fabs(r_disp(q)-L)){
                      r_disp_min(q) = r_disp(q);
                   }
                }
                if(fabs(r_disp(q)+L) < fabs(r_disp(q)-L)){
                   if(fabs(r_disp(q)+L) < fabs(r_disp(q))){
                      r_disp_min(q) = r_disp(q)+L;
                   }
                }
              }

            r_disp_2 = r_disp_min(0)*r_disp_min(0) + r_disp_min(1)*r_disp_min(1) + r_disp_min(2)*r_disp_min(2);

            r_sum = r_sum + r_disp_2;


//            for(j=0;j<dim;j++){
//                myfile << " " << r(j);
//            }
//            for(k=0;k<dim;k++){
//                myfile << " " << v(k);
//            }

//            myfile << "\n";


            v_half = v + F*dt/2;
            r_new = r + v_half*dt;

            for(j=0;j<dim;j++){  // periodic boundary conditions


                r_new(j) = fmod(r_new(j), (L));

                if (r_new(j)<0){
                    r_new(j) = L + r_new(j);
                }
             }

             atoms[i].setPosition(r_new);

        }

        E_tot = E_k + pot;

        T = 2*E_k/(3*numpart);

        P = numpart/V*T + 1/(3*V)*Fr_sum;

        mean_disp = 1.0/numpart*r_sum;
        D = mean_disp/(6*t*dt);

//        myfile << E_tot << "\n";
//        myfile << T << "\n";
//        myfile << P;
//        myfile << " " << T << "\n";
//        myfile << mean_disp << " " << D << " " << T << "\n";
//        myfile << E_tot << " " << T << "\n";

        if(t>200){
            for(i=0;i<numpart;i++){  // Finner g(r)
                r = atoms[i].getPosition();
                for(j=0;j<numpart;j++){
                    if(i!=j){
                        r2 = atoms[j].getPosition();
                        d = r - r2;

                        for(q=0;q<dim;q++){   // minimum image convention
                            if(fabs(d(q)-L) < fabs(d(q))){
                                if(fabs(d(q)-L) < fabs(d(q)+L)){
                                    r_min(q) = d(q)-L;
                                }
                            }
                            if(fabs(d(q)) < fabs(d(q)+L)){
                                if(fabs(d(q)) < fabs(d(q)-L)){
                                    r_min(q) = d(q);
                                }
                            }
                            if(fabs(d(q)+L) < fabs(d(q)-L)){
                                if(fabs(d(q)+L) < fabs(d(q))){
                                    r_min(q) = d(q)+L;
                                }
                            }
                        }

                        d_length = r_min(0)*r_min(0) + r_min(1)*r_min(1) + r_min(2)*r_min(2);

                        if(d_length <= L/2.0){
                            dist = int (d_length/(L/2.0)*bins);
                            volume = 4.0/3.0*pi*pow ((dist+1.0)*L/2.0/bins,3) - 4.0/3.0*pi*pow (dist*L/2.0/bins,3);
                            g_r(dist) = g_r(dist) + 1.0/volume;
                        }
                    }
                }
            }
        }


        for(i=0;i<bins;i++){
            myfile << g_r(i) << " ";
        }
        myfile << "\n";

//        for(i=0;i<numpart;i++){  // Berendsen
//            v = atoms[i].getVelocity();
//            gamma = sqrt(1 + dt/tau*(T_bath/T - 1));
//            for(j=0;j<dim;j++){
//                v(j) = v(j)*gamma;
//            }
//            atoms[i].setVelocity(v);
//        }

        if(t<200){
            for(i=0;i<numpart;i++){  // Andersen
                v = atoms[i].getVelocity();
                tall = randu();
                if(tall<(dt/tau)){
                    v_new = randn(3)*sqrt(T_bath);
                    atoms[i].setVelocity(v_new);
                }
            }
        }


        E_k = 0.0;
        pot = 0.0;
        Fr_sum = 0.0;
        r_sum = 0.0;

        Cell_container Con(cells_x);


        for(i=0;i<numpart;i++){  // Setter nye posisjoner i celler
            r = atoms[i].getPosition();
            nx = int (r(0)/r_cut);
            ny = int (r(1)/r_cut);
            nz = int (r(2)/r_cut);

            if(nx>cells_x-1){
                nx = cells_x-1;
            }
            if(ny>cells_x-1){
                ny = cells_x-1;
            }
            if(nz>cells_x-1){
                nz = cells_x-1;
            }


            Con.container(nx,ny,nz,&atoms[i]);

        }



        for(i=0;i<cells_x;i++){  // Finn kraft og fart
            for(j=0;j<cells_x;j++){
                for(k=0;k<cells_x;k++){
                    celle = Con.getContainer(i,j,k);

                    for(l=0;l<celle.size();l++){

                        r = celle[l]->getPosition();

                        for(g=0;g<dim;g++){
                            F(g) = 0;
                        }

                        for(g=0;g<celle.size();g++){   // egen celle
                            if(g!=l){
                                r2 = celle[g]->getPosition();

                                d = r - r2;
                                d_length = sqrt(d(0)*d(0) + d(1)*d(1) + d(2)*d(2));

                                pot = pot + 4*(pow (1/d_length,12) - pow (1/d_length,6));

                                for(h=0;h<dim;h++){
                                    F(h) = F(h) + 24*(2*pow (1/d_length,12) - pow (1/d_length,6))*d(h)/(d_length*d_length);
                                    Fr_sum = Fr_sum + F(h)*d(h);
                                }
                            }
                        }

                        for(g=-1;g<2;g++){   // naboceller
                            for(h=-1;h<2;h++){
                                for(o=-1;o<2;o++){
                                    if(!(g==0 && h==0 && o==0)){

                                        i2 = i+g;
                                        j2 = j+h;
                                        k2 = k+o;

                                        if((i+g) < 0){
                                            i2 = cells_x-1;
                                        }
                                        if((j+h) < 0){
                                            j2 = cells_x-1;
                                        }
                                        if((k+o) < 0){
                                            k2 = cells_x-1;
                                        }

                                        if((i+g) == cells_x){
                                            i2 = 0;
                                        }
                                        if((j+h) == cells_x){
                                            j2 = 0;
                                        }
                                        if((k+o) == cells_x){
                                            k2 = 0;
                                        }

                                        celle2 = Con.getContainer(i2,j2,k2);


                                        for(p=0;p<celle2.size();p++){
                                            r2 = celle2[p]->getPosition();
                                            d = r-r2;

                                            for(q=0;q<dim;q++){    // minimum image convention
                                                if(fabs(d(q)-L) < fabs(d(q))){
                                                   if(fabs(d(q)-L) < fabs(d(q)+L)){
                                                       r_min(q) = d(q)-L;
                                                   }
                                                }
                                                if(fabs(d(q)) < fabs(d(q)+L)){
                                                   if(fabs(d(q)) < fabs(d(q)-L)){
                                                      r_min(q) = d(q);
                                                   }
                                                }
                                                if(fabs(d(q)+L) < fabs(d(q)-L)){
                                                   if(fabs(d(q)+L) < fabs(d(q))){
                                                      r_min(q) = d(q)+L;
                                                   }
                                                }
                                            }


                                            d_length = sqrt(r_min(0)*r_min(0) + r_min(1)*r_min(1) + r_min(2)*r_min(2));

                                            pot = pot + 4*(pow (1/d_length,12) - pow (1/d_length,6));

                                            for(q=0;q<dim;q++){
                                                F(q) = F(q) + 24*(2*pow (1/d_length,12) - pow (1/d_length,6))*r_min(q)/(d_length*d_length);
                                                Fr_sum = Fr_sum + F(q)*r_min(q);
                                            }
                                        }
                                    }
                                }
                            }
                        }


                        v = celle[l]->getVelocity();
                        F_old = celle[l]->getForce();
                        v_half = v + F_old*dt/2;
                        v_new = v_half + F*dt/2;

                        celle[l]->setForce(F);
                        celle[l]->setVelocity(v_new);

                    }
                }
            }
        }


//        myfile.close();

    }
    myfile.close();
}







//for(i=0;i<numpart;i++){  // Finn kraft første gang

//    r = atoms[i].getPosition();

//    for(j=0;j<dim;j++){
//        F(j) = 0;
//    }

//    for(j=0;j<numpart;j++){
//        if (i!=j){
//            r2 = atoms[j].getPosition();
//            d = r - r2;

//            for(k=0;k<dim;k++){
//                if(abs(d(k)-L) < abs(d(k))){
//                    if(abs(d(k)-L) < abs(d(k)+L)){
//                        r_min(k) = d(k)-L;
//                    }
//                }
//                if(abs(d(k)) < abs(d(k)+L)){
//                    if(abs(d(k)) < abs(d(k)-L)){
//                        r_min(k) = d(k);
//                    }
//                }
//                if(abs(d(k)+L) < abs(d(k)-L)){
//                    if(abs(d(k)+L) < abs(d(k))){
//                        r_min(k) = d(k)+L;
//                    }
//                }
//            }

//            d_length = sqrt(r_min(0)*r_min(0) + r_min(1)*r_min(1) + r_min(2)*r_min(2));

//            for(k=0;k<dim;k++){
//                F(k) = F(k) + 24*(2*pow (1/d_length,12) - pow (1/d_length,6))*r_min(k)/(d_length*d_length);
//            }
//        }

//    }

//    atoms[i].setForce(F);


//}



//for(j=0;j<numpart;j++){  // Finn kraft  og fart

//    r = atoms[j].getPosition();

//    for(k=0;k<dim;k++){
//            F(k) = 0;
//    }

//    for(l=0;l<numpart;l++){
//        if (l!=j){
//            r2 = atoms[l].getPosition();

//            d = r - r2;

//            for(k=0;k<dim;k++){
//                if(abs(d(k)-L) < abs(d(k))){
//                    if(abs(d(k)-L) < abs(d(k)+L)){
//                        r_min(k) = d(k)-L;
//                    }
//                }
//                if(abs(d(k)) < abs(d(k)+L)){
//                    if(abs(d(k)) < abs(d(k)-L)){
//                        r_min(k) = d(k);
//                    }
//                }
//                if(abs(d(k)+L) < abs(d(k)-L)){
//                    if(abs(d(k)+L) < abs(d(k))){
//                        r_min(k) = d(k)+L;
//                    }
//                }
//            }



//            d_length = sqrt(r_min(0)*r_min(0) + r_min(1)*r_min(1) + r_min(2)*r_min(2));


//            for(k=0;k<dim;k++){
//                F(k) = F(k) + 24*(2*pow (1/d_length,12) - pow (1/d_length,6))*r_min(k)/(d_length*d_length);
//                //F[k] = F[k] - 24*eps*(2*pow (sig/d_length,12) - pow (sig/d_length,6))*d[k]/(d_length*d_length);
//            }

//        }

//    }

//    v = atoms[j].getVelocity();
//    F_old = atoms[j].getForce();
//    v_half = v + F_old*dt/2;
//    v_new = v_half + F*dt/2;

//    atoms[j].setForce(F);
//    atoms[j].setVelocity(v_new);
//    //cout << F << endl;

//}
//}
//myfile.close();
//}
