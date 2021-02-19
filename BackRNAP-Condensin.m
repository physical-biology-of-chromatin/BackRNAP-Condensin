% Gillespie simulations of the backtrack model of Rivosecchi et al
% developed by Daniel Jost (CNRS, contact: daniel.jost@ens-lyon.fr)
% 2021 - 02 - 15

% the simple model may be obtained with term=termb >0 and koffmin=koffmax=0

%% head to tail situation
clear all
N=100;%size of the region 1 bin=100bp
Nterm=20;%start of the termination zone
init=0.05;%RNAP initiation rate (1/s)
v=0.4; %speed rate of mobile RNAP in the gene body (in bin/s)
vterm=v; %speed rate of mobile RNAP in the termination zone (in bin/s)
term=0.; %termination rate for mobile RNAP
termb=0.2; %termination rate for backtracked RNAP

kon=0.2; %rate to switch from backtracked to mobile state (1/s)
koffmin=0.002;  %minimal rate to switch from mobile to backtracked state in the gene body (1/s)
koffmax=200*koffmin; %maximal rate to switch from mobile to backtracked state in the termination zone (1/s)
koff=koffmin*ones(1,N); % rate to switch from mobile to backtracked state along the region
koff([Nterm:(Nterm+3)])=koffmin+(koffmax-koffmin)*([1:4])/4;
koff([(Nterm+4):N])=koffmax;

vc=10;%speed rate of condensin (in bin/s)
jump=0.15; %speed rate to jump over mobile RNAP (in bin/s)
jumpb=0.; %speed rate to jump over backtracked RNAP (in bin/s)

ntraj=100000; % number of simulated trajectories
Teq=Nterm/v*10; %time to equilibrate the RNAP organization (in s) 
statetot=zeros(ntraj,N); % output 1 : RNAP organization after equilibration; statetot(i,:) represents the RNAP occupancy of trajectory i at time Teq along the region; statetot(i,j)=0 if no RNAP at bin j, =1 if mobile RNAP, =2 if backtraked RNAP.
dt=zeros(ntraj,N); % output 2: residence time of conensin; dt(i,j) represents the time spent by condensin at bin j during its travel across the region simulated in trajectory i.

tic
for it=1:ntraj
    %convergence to steady-state for the tasep model
    t=0;
    state=zeros(1,N);
    ri=init;
    rv=zeros(1,N);
    while t<Teq
        s=state;
        ro=rv;
        rtot=ri+sum(rv);
        t=t-log(rand(1))/rtot;
        r=rand(1);
        if r<=ri/rtot
            ri=0;
            state(1)=1;rv(1)=koff(1);
            if state(2)==0
                rv(1)=rv(1)+v;
            end      
        else
            a=find(r<=(ri+cumsum(rv))/rtot);
            l=a(1);
            if l==N 
                state(N)=0;
                rv(N)=0;
                if state(N-1)==1
                    rv(N-1)=rv(N-1)+vterm;
                end
            else
                if (state(l)==2) 
                    if l>=Nterm 
                        if rand(1)<termb/rv(l) 
                            state(l)=0;
                            rv(l)=0;
                            if state(l-1)==1
                                if (l-1)>=Nterm
                                    rv(l-1)=rv(l-1)+vterm;
                                else
                                    rv(l-1)=rv(l-1)+v;
                                end
                            end
                        else 
                            state(l)=1;
                            rv(l)=koff(l)+term;
                            if state(l+1)==0
                                rv(l)=rv(l)+vterm;
                            end
                        end
                    else               
                        state(l)=1;
                        rv(l)=koff(l);
                        if state(l+1)==0
                           rv(l)=rv(l)+v;
                        end
                    end
                    
                else 
                    e=0;r1=rand(1);
                    if l>=Nterm 
                       if (r1<term/rv(l)) 
                          state(l)=0;
                          rv(l)=0;
                          if state(l-1)==1
                             if (l-1)>=Nterm
                                 rv(l-1)=rv(l-1)+vterm;
                             else
                                 rv(l-1)=rv(l-1)+v;
                             end
                          end
                       else
                           e=1;r1=r1-term/rv(l);
                       end
                    end
                    if (l<Nterm)||(e==1)  
                        if r1<koff(l)/rv(l) 
                            state(l)=2;
                            rv(l)=kon;
                            if (l>=Nterm)
                                rv(l)=rv(l)+termb;
                            end
                        else 
                            state(l+1)=1;
                            state(l)=0;
                            rv(l)=0;
                            rv(l+1)=koff(l+1);
                            if (l+1)>=Nterm
                                rv(l+1)=rv(l+1)+term;
                            end
                            if (l+2<=N)&& state(l+2)==0
                                if l+1>=Nterm
                                    rv(l+1)=rv(l+1)+vterm;
                                else
                                    rv(l+1)=rv(l+1)+v;
                                end
                            end
                            
                            if l==1
                                ri=init;
                            elseif state(l-1)==1
                                if (l-1)>=Nterm
                                    rv(l-1)=rv(l-1)+vterm;
                                else
                                    rv(l-1)=rv(l-1)+v;
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    statetot(it,:)=s;
    t=0;to=0;
    state=s;
    rv=ro;
    if state(1)==0
        ri=init;
        rc=vc;
    else
        ri=0;
        if state(1)==1
            rc=jump;
        elseif state(1)==2
            rc=jumpb;
        end
    end
    
    p=1;
    while p<=N
        rtot=rc+ri+sum(rv);
        t=t-log(rand(1))/rtot;
        r=rand(1);
        if r<=rc/rtot
            dt(it,p)=dt(it,p)+t-to;to=t;
            p=p+1;
            if (p<=N)
                if state(p)==0
                    rc=vc;
                elseif state(p)==1
                    rc=jump;
                elseif state(p)==2
                    rc=jumpb;
                end
            end
        elseif r<=(ri+rc)/rtot
            ri=0;
            state(1)=1;rv(1)=koff(1);
            if state(2)==0
                rv(1)=rv(1)+v;
            end
            if p==1
                rc=jump;
            end
        else
            a=find(r<=(ri+rc+cumsum(rv))/rtot);
            l=a(1);
            if l==N  
                state(N)=0;
                rv(N)=0;
                if state(N-1)==1
                    rv(N-1)=rv(N-1)+vterm;
                end
                if p==N
                    rc=vc;
                end
            else
                if (state(l)==2) 
                    if l>=Nterm 
                        if rand(1)<termb/rv(l) 
                            state(l)=0;
                            rv(l)=0;
                            if state(l-1)==1
                                if (l-1)>=Nterm
                                    rv(l-1)=rv(l-1)+vterm;
                                else
                                    rv(l-1)=rv(l-1)+v;
                                end
                            end
                            if p==l
                                rc=vc;
                            end
                        else 
                            state(l)=1;
                            rv(l)=koff(l)+term;
                            if state(l+1)==0
                                rv(l)=rv(l)+vterm;
                            end
                            if p==l
                                rc=jump;
                            end
                        end
                    else               
                        state(l)=1;
                        rv(l)=koff(l);
                        if state(l+1)==0
                           rv(l)=rv(l)+v;
                        end
                        if p==l
                            rc=jump;
                        end
                    end
                    
                else 
                    e=0;r1=rand(1);
                    if l>=Nterm 
                       if (r1<term/rv(l)) 
                          state(l)=0;
                          rv(l)=0;
                          if state(l-1)==1
                             if (l-1)>=Nterm
                                 rv(l-1)=rv(l-1)+vterm;
                             else
                                 rv(l-1)=rv(l-1)+v;
                             end
                          end
                          if p==l
                              rc=vc;
                          end
                       else
                           e=1;r1=r1-term/rv(l);
                       end
                    end
                    if (l<Nterm)||(e==1) 
                        if r1<koff(l)/rv(l) 
                            state(l)=2;
                            rv(l)=kon;
                            if (l>=Nterm)
                                rv(l)=rv(l)+termb;
                            end
                            if (p==l)
                                rc=jumpb;
                            end
                        else 
                            state(l+1)=1;
                            state(l)=0;
                            rv(l)=0;
                            rv(l+1)=koff(l+1);
                            if (l+1)>=Nterm
                                rv(l+1)=rv(l+1)+term;
                            end
                            if (l+2<=N)&& state(l+2)==0
                                if l+1>=Nterm
                                    rv(l+1)=rv(l+1)+vterm;
                                else
                                    rv(l+1)=rv(l+1)+v;
                                end
                            end
                            
                            if l==1
                                ri=init;
                            elseif state(l-1)==1
                                if (l-1)>=Nterm
                                    rv(l-1)=rv(l-1)+vterm;
                                else
                                    rv(l-1)=rv(l-1)+v;
                                end
                            end
                            if p==l
                                rc=vc;
                            elseif p==l+1
                                dt(it,p)=dt(it,p)+t-to;to=t;
                                p=p+1;
                                if (p<=N)
                                    if state(p)==0
                                        rc=vc;
                                    elseif state(p)==1
                                        rc=jump;
                                    elseif state(p)==2
                                        rc=jumpb;
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    
end
toc

%% head to head situation
clear all
N=100;%size of the region 1 bin=100bp
Nterm=20;%start of the termination zone
init=0.05;%RNAP initiation rate (1/s)
v=0.4; %speed rate of mobile RNAP in the gene body (in bin/s)
vterm=v; %speed rate of mobile RNAP in the termination zone (in bin/s)
term=0.; %termination rate for mobile RNAP
termb=0.2; %termination rate for backtracked RNAP

kon=0.2; %rate to switch from backtracked to mobile state (1/s)
koffmin=0.002;  %minimal rate to switch from mobile to backtracked state in the gene body (1/s)
koffmax=200*koffmin; %maximal rate to switch from mobile to backtracked state in the termination zone (1/s)
koff=koffmin*ones(1,N); % rate to switch from mobile to backtracked state along the region
koff([Nterm:(Nterm+3)])=koffmin+(koffmax-koffmin)*([1:4])/4;
koff([(Nterm+4):N])=koffmax;

vc=10;%speed rate of condensin (in bin/s)
jump=0.15; %speed rate to jump over mobile RNAP (in bin/s)
jumpb=0.; %speed rate to jump over backtracked RNAP (in bin/s)

ntraj=100000; % number of simulated trajectories
Teq=Nterm/v*10; %time to equilibrate the RNAP organization (in s) 
statetot=zeros(ntraj,N); % output 1 : RNAP organization after equilibration; statetot(i,:) represents the RNAP occupancy of trajectory i at time Teq along the region; statetot(i,j)=0 if no RNAP at bin j, =1 if mobile RNAP, =2 if backtraked RNAP.
dt=zeros(ntraj,N); % output 2: residence time of conensin; dt(i,j) represents the time spent by condensin at bin j during its travel across the region simulated in trajectory i.

tic
for it=1:ntraj
    t=0;
    state=zeros(1,N);
    ri=init;
    rv=zeros(1,N);
    while t<Teq
        s=state;
        ro=rv;
        rtot=ri+sum(rv);
        t=t-log(rand(1))/rtot;
        r=rand(1);
        if r<=ri/rtot
            ri=0;
            state(1)=1;rv(1)=koff(1);
            if state(2)==0
                rv(1)=rv(1)+v;
            end      
        else
            a=find(r<=(ri+cumsum(rv))/rtot);
            l=a(1);
            if l==N  
                state(N)=0;
                rv(N)=0;
                if state(N-1)==1
                    rv(N-1)=rv(N-1)+vterm;
                end
            else
                if (state(l)==2) 
                    if l>=Nterm 
                        if rand(1)<termb/rv(l)
                            state(l)=0;
                            rv(l)=0;
                            if state(l-1)==1
                                if (l-1)>=Nterm
                                    rv(l-1)=rv(l-1)+vterm;
                                else
                                    rv(l-1)=rv(l-1)+v;
                                end
                            end
                        else 
                            state(l)=1;
                            rv(l)=koff(l)+term;
                            if state(l+1)==0
                                rv(l)=rv(l)+vterm;
                            end
                        end
                    else                  
                        state(l)=1;
                        rv(l)=koff(l);
                        if state(l+1)==0
                           rv(l)=rv(l)+v;
                        end
                    end
                    
                else 
                    e=0;r1=rand(1);
                    if l>=Nterm 
                       if (r1<term/rv(l)) 
                          state(l)=0;
                          rv(l)=0;
                          if state(l-1)==1
                             if (l-1)>=Nterm
                                 rv(l-1)=rv(l-1)+vterm;
                             else
                                 rv(l-1)=rv(l-1)+v;
                             end
                          end
                       else
                           e=1;r1=r1-term/rv(l);
                       end
                    end
                    if (l<Nterm)||(e==1)  
                        if r1<koff(l)/rv(l) 
                            state(l)=2;
                            rv(l)=kon;
                            if (l>=Nterm)
                                rv(l)=rv(l)+termb;
                            end
                        else 
                            state(l+1)=1;
                            state(l)=0;
                            rv(l)=0;
                            rv(l+1)=koff(l+1);
                            if (l+1)>=Nterm
                                rv(l+1)=rv(l+1)+term;
                            end
                            if (l+2<=N)&& state(l+2)==0
                                if l+1>=Nterm
                                    rv(l+1)=rv(l+1)+vterm;
                                else
                                    rv(l+1)=rv(l+1)+v;
                                end
                            end
                            
                            if l==1
                                ri=init;
                            elseif state(l-1)==1
                                if (l-1)>=Nterm
                                    rv(l-1)=rv(l-1)+vterm;
                                else
                                    rv(l-1)=rv(l-1)+v;
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    statetot(it,:)=s;
    t=0;to=0;
    state=s;
    rv=ro;
    if state(1)==0
        ri=init;
    else
        ri=0;
    end
    
    p=N+1;
    if state(N)==0
        rc=vc;
    elseif state(N)==1
        rc=jump;
    elseif state(N)==2
        rc=jumpb;
    end
    
    while p>=2
        rtot=rc+ri+sum(rv);
        t=t-log(rand(1))/rtot;
        r=rand(1);
        if r<=rc/rtot
            dt(it,p-1)=dt(it,p-1)+t-to;to=t;
            p=p-1;
            if (p>=2)
                if state(p-1)==0
                    rc=vc;
                elseif state(p-1)==1
                    rc=jump;
                elseif state(p-1)==2
                    rc=jumpb;
                end
            end
        elseif r<=(ri+rc)/rtot
            ri=0;
            state(1)=1;rv(1)=koff(1);
            if state(2)==0
                rv(1)=rv(1)+v;
            end
            if p==2
                rc=jump;
            end
        else
            a=find(r<=(ri+rc+cumsum(rv))/rtot);
            l=a(1);
            if l==N  
                state(N)=0;
                rv(N)=0;
                if state(N-1)==1
                    rv(N-1)=rv(N-1)+vterm;
                end
                if p==N+1
                    rc=vc;
                end
            else
                if (state(l)==2) 
                    if l>=Nterm 
                        if rand(1)<termb/rv(l) 
                            state(l)=0;
                            rv(l)=0;
                            if state(l-1)==1
                                if (l-1)>=Nterm
                                    rv(l-1)=rv(l-1)+vterm;
                                else
                                    rv(l-1)=rv(l-1)+v;
                                end
                            end
                            if p==l+1
                                rc=vc;
                            end
                        else 
                            state(l)=1;
                            rv(l)=koff(l)+term;
                            if state(l+1)==0
                                rv(l)=rv(l)+vterm;
                            end
                            if p==l+1
                                rc=jump;
                            end
                        end
                    else               
                        state(l)=1;
                        rv(l)=koff(l);
                        if state(l+1)==0
                           rv(l)=rv(l)+v;
                        end
                        if p==l+1
                            rc=jump;
                        end
                    end
                    
                else 
                    e=0;r1=rand(1);
                    if l>=Nterm 
                       if (r1<term/rv(l))
                          state(l)=0;
                          rv(l)=0;
                          if state(l-1)==1
                             if (l-1)>=Nterm
                                 rv(l-1)=rv(l-1)+vterm;
                             else
                                 rv(l-1)=rv(l-1)+v;
                             end
                          end
                          if p==l+1
                              rc=vc;
                          end
                       else
                           e=1;r1=r1-term/rv(l);
                       end
                    end
                    if (l<Nterm)||(e==1)  
                        if r1<koff(l)/rv(l)
                            state(l)=2;
                            rv(l)=kon;
                            if (l>=Nterm)
                                rv(l)=rv(l)+termb;
                            end
                            if (p==l+1)
                                rc=jumpb;
                            end
                        else 
                            state(l+1)=1;
                            state(l)=0;
                            rv(l)=0;
                            rv(l+1)=koff(l+1);
                            if (l+1)>=Nterm
                                rv(l+1)=rv(l+1)+term;
                            end
                            if (l+2<=N)&& state(l+2)==0
                                if l+1>=Nterm
                                    rv(l+1)=rv(l+1)+vterm;
                                else
                                    rv(l+1)=rv(l+1)+v;
                                end
                            end
                            
                            if l==1
                                ri=init;
                            elseif state(l-1)==1
                                if (l-1)>=Nterm
                                    rv(l-1)=rv(l-1)+vterm;
                                else
                                    rv(l-1)=rv(l-1)+v;
                                end
                            end
                            if p==l+1
                                dt(it,p-1)=dt(it,p-1)+t-to;to=t;
                                p=p+1;
                                rc=jump;
                            elseif (l<=N-1) && p==l+2
                                rc=jump;
                            end
                        end
                    end
                end
            end
        end
    end
    
end
toc