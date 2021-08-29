#!/usr/bin/env python
# coding: utf-8

# In[ ]:


#3D Cavitry Staggered grids (uniform grids)
import numpy as np
from timeit import default_timer as timer
#import cupy as cp

#Global method
def TDMA(a,b,c,d,g,factor,nx,ny,nz,x,y,z):
    if x == 1: #x-direction TDMA
        for i in range(1,nx): #remove of a coefficient
            
            factor[i,:,:,:] = a[i,:,:,:]/b[i-1,:,:,:] #factor,a는 0 index는 버리는 것
            b[i,:,:,:] = b[i,:,:,:] - factor[i,:,:,:]*c[i-1,:,:,:]
            d[i,:,:,:] = d[i,:,:,:] - factor[i,:,:,:]*d[i-1,:,:,:]
            
        g[nx-1,:,:,:] = d[nx-1,:,:,:]/b[nx-1,:,:,:]
        
        for i in range(nx-2,(-1),-1):
            g[i,:,:,:] = (d[i,:,:,:]-c[i,:,:,:]*g[i+1,:,:,:])/b[i,:,:,:]
            
    elif y == 1: #y-direction TDMA
        for j in range(1,ny): #remove of a coefficient
            
            factor[:,j,:,:] = a[:,j,:,:]/b[:,j-1,:,:] 
            b[:,j,:,:] = b[:,j,:,:] - factor[:,j,:,:]*c[:,j-1,:,:]
            d[:,j,:,:] = d[:,j,:,:] - factor[:,j,:,:]*d[:,j-1,:,:]
            
        g[:,ny-1,:,:] = d[:,ny-1,:,:]/b[:,ny-1,:,:]
        
        for j in range(ny-2,(-1),-1):
            g[:,j,:,:] = (d[:,j,:,:]-c[:,j,:,:]*g[:,j+1,:,:])/b[:,j,:,:]
    
    elif z == 1: #y-direction TDMA
        for k in range(1,nz): #remove of a coefficient
            
            factor[:,:,k,:] = a[:,:,k,:]/b[:,:,k-1,:] 
            b[:,:,k,:] = b[:,:,k,:] - factor[:,:,k,:]*c[:,:,k-1,:]
            d[:,:,k,:] = d[:,:,k,:] - factor[:,:,k,:]*d[:,:,k-1,:]
            
        g[:,:,nz-1,:] = d[:,:,nz-1,:]/b[:,:,nz-1,:]
        
        for k in range(nz-2,(-1),-1):
            g[:,:,k,:] = (d[:,:,k,:]-c[:,:,k,:]*g[:,:,k+1,:])/b[:,:,k,:]
    
    return g

def Poisson(p_n,p_o,u_mid,v_mid,w_mid,dt,nx,ny,nz,dx,dy,dz,epsil_p):
    
    steps = int(0)
    
    while(1):
        
        sum_error = np.zeros([1],dtype = np.float64)
        
        if steps ==0:
            p_o[0,1:ny,1:nz,:] = p_o[1,1:ny,1:nz,:]; p_o[nx,1:ny,1:nz,:] = p_o[nx-1,1:ny,1:nz,:]
            p_o[1:nx,0,1:nz,:] = p_o[1:nx,1,1:nz,:]; p_o[1:nx,ny,1:nz,:] = p_o[1:nx,ny-1,1:nz,:]  
            p_o[1:nx,1:ny,0,:] = p_o[1:nx,1:ny,1,:]; p_o[1:nx,1:ny,nz,:] = p_o[1:nx,1:ny,nz-1,:]
        
        p_n[1:nx,1:ny,1:nz,:] = (p_o[2:nx+1,1:ny,1:nz,:] + p_o[0:nx-1,1:ny,1:nz,:]+ p_o[1:nx,2:ny+1,1:nz,:] + p_o[1:nx,0:ny-1,1:nz,:] +                               p_o[1:nx,1:ny,2:nz+1,:] + p_o[1:nx,1:ny,0:nz-1,:]                               - dx**2/dt*((u_mid[2:nx+1,1:ny,1:nz]-u_mid[1:nx,1:ny,1:nz])/dx                               +(v_mid[1:nx,2:ny+1,1:nz]-v_mid[1:nx,1:ny,1:nz])/dy + (w_mid[1:nx,1:ny,2:nz+1]-w_mid[1:nx,1:ny,1:nz])/dz))/6.0
        
        #Neumman B.C
        p_n[0,1:ny,1:nz,:] = p_n[1,1:ny,1:nz,:]; p_n[nx,1:ny,1:nz,:] = p_n[nx-1,1:ny,1:nz,:]
        p_n[1:nx,0,1:nz,:] = p_n[1:nx,1,1:nz,:]; p_n[1:nx,ny,1:nz,:] = p_n[1:nx,ny-1,1:nz,:]
        p_n[1:nx,1:ny,0,:] = p_n[1:nx,1:ny,1,:]; p_n[1:nx,1:ny,nz,:] = p_n[1:nx,1:ny,nz-1,:]
        
        sum_error[:] = np.sum(abs(p_n[1:nx,1:ny,1:nz] - p_o[1:nx,1:ny,1:nz]))/((nx-1)*(ny-1)*(nz-1))
        
        p_o[:]= p_n[:]
        
        #print(sum_error)
        if (sum_error<epsil_p):
            break
            
        steps = steps+1
        
    return p_o, steps

def main():
    
    #grids shape
    nx = 64; ny = 64; nz = 64
    
    #---------------------- Define variables ---------------------#
    #Dfine 3-D velocity
    u_o = np.zeros([nx+1,ny+1,nz+1,1], dtype = np.float64)
    u_n = np.zeros([nx+1,ny+1,nz+1,1], dtype = np.float64)    
    v_o = np.zeros([nx+1,ny+1,nz+1,1], dtype = np.float64)
    v_n = np.zeros([nx+1,ny+1,nz+1,1], dtype = np.float64)
    w_o = np.zeros([nx+1,ny+1,nz+1,1], dtype = np.float64)
    w_n = np.zeros([nx+1,ny+1,nz+1,1], dtype = np.float64)
    
    u_c = np.zeros([nx,ny,nz,3], dtype = np.float64)
    stream = np.zeros([nx,ny,nz,1], dtype = np.float64)
    vorticity = np.zeros([nx-1,ny-1,nz-1,3],dtype = np.float64)
    p_o = np.zeros([nx+1,ny+1,nz+1,1], dtype = np.float64)
    p_n = np.zeros([nx+1,ny+1,nz+1,1], dtype = np.float64)
    
    #Coefficient
    #x-direction coefficient
    a_x = np.zeros([nx+1,1,1,1],dtype = np.float64); b_x = np.zeros([nx+1,1,1,1],dtype = np.float64)
    c_x = np.zeros([nx+1,1,1,1],dtype = np.float64); d_u = np.zeros([nx-2,ny-1,nz-1,1],dtype = np.float64)
    factor_x = np.zeros([nx+1,1,1,1],dtype = np.float64)

    #y-direction coefficient
    a_y = np.zeros([1,ny+1,1,1],dtype = np.float64); b_y = np.zeros([1,ny+1,1,1],dtype = np.float64)
    c_y = np.zeros([1,ny+1,1,1],dtype = np.float64); d_v = np.zeros([nx-1,ny-2,nz-1,1],dtype = np.float64)
    factor_y = np.zeros([1,ny+1,1,1],dtype = np.float64)

    #z-direction coefficient
    a_z = np.zeros([1,1,nz+1,1],dtype = np.float64); b_z = np.zeros([1,1,nz+1,1],dtype = np.float64)
    c_z = np.zeros([1,1,nz+1,1],dtype = np.float64); d_w = np.zeros([nx-1,ny-1,nz-2,1],dtype = np.float64)
    factor_z = np.zeros([1,1,nz+1,1],dtype = np.float64)
    
    #intermediate factor for TDMA
    g = np.zeros([nx+1,ny+1,nz+1,1], dtype = np.float64)
    h = np.zeros([nx+1,ny+1,nz+1,1], dtype = np.float64)
    q = np.zeros([nx+1,ny+1,nz+1,1], dtype = np.float64)
    u_mid = np.zeros([nx+1,ny+1,nz+1,1], dtype = np.float64)
    v_mid = np.zeros([nx+1,ny+1,nz+1,1], dtype = np.float64)
    w_mid = np.zeros([nx+1,ny+1,nz+1,1], dtype = np.float64)
    
    #gradient of velocity
    dudx = np.zeros([nx+1,ny+1,nz+1,1], dtype = np.float64)
    dudy = np.zeros([nx+1,ny+1,nz+1,1], dtype = np.float64)    
    dudz = np.zeros([nx+1,ny+1,nz+1,1], dtype = np.float64)  
    dvdx = np.zeros([nx+1,ny+1,nz+1,1], dtype = np.float64)
    dvdy = np.zeros([nx+1,ny+1,nz+1,1], dtype = np.float64)
    dvdz = np.zeros([nx+1,ny+1,nz+1,1], dtype = np.float64)
    dwdx = np.zeros([nx+1,ny+1,nz+1,1], dtype = np.float64)
    dwdy = np.zeros([nx+1,ny+1,nz+1,1], dtype = np.float64)
    dwdz = np.zeros([nx+1,ny+1,nz+1,1], dtype = np.float64)
    
    dudx_2 = np.zeros([nx+1,ny+1,nz+1,1], dtype = np.float64)
    dudy_2 = np.zeros([nx+1,ny+1,nz+1,1], dtype = np.float64) 
    dudz_2 = np.zeros([nx+1,ny+1,nz+1,1], dtype = np.float64) 
    dvdx_2 = np.zeros([nx+1,ny+1,nz+1,1], dtype = np.float64)
    dvdy_2 = np.zeros([nx+1,ny+1,nz+1,1], dtype = np.float64)
    dvdz_2 = np.zeros([nx+1,ny+1,nz+1,1], dtype = np.float64)
    dwdx_2 = np.zeros([nx+1,ny+1,nz+1,1], dtype = np.float64)
    dwdy_2 = np.zeros([nx+1,ny+1,nz+1,1], dtype = np.float64)
    dwdz_2 = np.zeros([nx+1,ny+1,nz+1,1], dtype = np.float64)
    
    #Non-linear
    H_uo = np.zeros([nx+1,ny+1,nz+1,1], dtype = np.float64)
    H_vo = np.zeros([nx+1,ny+1,nz+1,1], dtype = np.float64)
    H_wo = np.zeros([nx+1,ny+1,nz+1,1], dtype = np.float64)
    H_up = np.zeros([nx+1,ny+1,nz+1,1], dtype = np.float64)
    H_vp = np.zeros([nx+1,ny+1,nz+1,1], dtype = np.float64)
    H_wp = np.zeros([nx+1,ny+1,nz+1,1], dtype = np.float64)
    #-------------------------------------------------------------#
    
    # Simulation parameter
    Re = 1000
    dt = 0.0025
    t_real =  int(0)
    epsil_p = 10**(-8)
    epsil_u = 10**(-8)
    
    #simulation factor
    x_len = 1.0; y_len = 1.0; z_len = 1.0
    
    dx = x_len/(nx-1); dy = y_len/(ny-1); dz = z_len/(nz-1) #delta panel 칸의 갯수로 나눠줘야함
    
    mu = dt/(2*Re*dx**2)
    
    start = timer()
    
    while(1):
        
        #update non-linear term
        H_up[2:nx,1:ny,1:nz] = H_uo[2:nx,1:ny,1:nz]
        H_vp[1:nx,2:ny,1:nz] = H_vo[1:nx,2:ny,1:nz]
        H_wp[1:nx,1:ny,2:nz] = H_wo[1:nx,1:ny,2:nz]
        
        #Boundary condition (자기 자신의 direction 외에는 모두 해줘야함)
        u_o[1:nx,ny,1:nz,:] = (2 - u_o[1:nx,ny-1,1:nz,:]); u_o[1:nx,0,1:nz,:] = - u_o[1:nx,1,1:nz,:]
        u_o[1:nx,1:ny,nz,:] =  -u_o[1:nx,1:ny,nz-1,:]; u_o[1:nx,1:ny,0,:] = - u_o[1:nx,1:ny,1,:]
        
        v_o[nx,1:ny,1:nz,:] = - v_o[nx-1,1:ny,1:nz,:]; v_o[0,1:ny,1:nz,:] = - v_o[1,1:ny,1:nz,:]
        v_o[1:nx,1:ny,nz,:] = - v_o[1:nx,1:ny,nz-1,:]; v_o[1:nx,1:ny,0,:] = - v_o[1:nx,1:ny,1,:]

        w_o[nx,1:ny,1:nz,:] = - w_o[nx-1,1:ny,1:nz,:]; w_o[0,1:ny,1:nz,:] = - w_o[1,1:ny,1:nz,:]
        w_o[1:nx,ny,1:nz,:] = - w_o[1:nx,ny-1,1:nz,:]; w_o[1:nx,0,1:nz,:] = - w_o[1:nx,1,1:nz,:]
    
        #Gradient
        dudx[2:nx,1:ny,1:nz,:] = 0.5*((u_o[3:nx+1,1:ny,1:nz,:] - u_o[2:nx,1:ny,1:nz,:])/dx +                                       (u_o[2:nx,1:ny,1:nz,:] - u_o[1:nx-1,1:ny,1:nz,:])/dx)
        dudy[1:nx,1:ny,1:nz,:] = 0.5*((u_o[1:nx,2:ny+1,1:nz,:] - u_o[1:nx,1:ny,1:nz,:])/dy +                                       (u_o[1:nx,1:ny,1:nz,:] - u_o[1:nx,0:ny-1,1:nz,:])/dy)
        dudz[1:nx,1:ny,1:nz,:] = 0.5*((u_o[1:nx,1:ny,2:nz+1,:] - u_o[1:nx,1:ny,1:nz,:])/dz +                                       (u_o[1:nx,1:ny,1:nz,:] - u_o[1:nx,1:ny,0:nz-1,:])/dz)

        dvdx[1:nx,1:ny,1:nz,:] = 0.5*((v_o[2:nx+1,1:ny,1:nz,:] - v_o[1:nx,1:ny,1:nz,:])/dx +                                       (v_o[1:nx,1:ny,1:nz,:] - v_o[0:nx-1,1:ny,1:nz,:])/dx)
        dvdy[1:nx,2:ny,1:nz,:] = 0.5*((v_o[1:nx,3:ny+1,1:nz,:] - v_o[1:nx,2:ny,1:nz,:])/dy +                                       (v_o[1:nx,2:ny,1:nz,:] - v_o[1:nx,1:ny-1,1:nz,:])/dy)
        dvdz[1:nx,1:ny,1:nz,:] = 0.5*((v_o[1:nx,1:ny,2:nz+1,:] - v_o[1:nx,1:ny,1:nz,:])/dz +                                       (v_o[1:nx,1:ny,1:nz,:] - v_o[1:nx,1:ny,0:nz-1,:])/dz)
        
        dwdx[1:nx,1:ny,1:nz,:] = 0.5*((w_o[2:nx+1,1:ny,1:nz,:] - w_o[1:nx,1:ny,1:nz,:])/dx +                                       (w_o[1:nx,1:ny,1:nz,:] - w_o[0:nx-1,1:ny,1:nz,:])/dx)
        dwdy[1:nx,1:ny,1:nz,:] = 0.5*((w_o[1:nx,2:ny+1,1:nz,:] - w_o[1:nx,1:ny,1:nz,:])/dy +                                       (w_o[1:nx,1:ny,1:nz,:] - w_o[1:nx,0:ny-1,1:nz,:])/dy)
        dwdz[1:nx,1:ny,2:nz,:] = 0.5*((w_o[1:nx,1:ny,3:nz+1,:] - w_o[1:nx,1:ny,2:nz,:])/dz +                                       (w_o[1:nx,1:ny,2:nz,:] - w_o[1:nx,1:ny,1:nz-1,:])/dz)        
        #dvdz, dwdz가 그냥 copy해서 dvdx, dwdx로 되어 있었음
        
        H_uo[2:nx,1:ny,1:nz] = u_o[2:nx,1:ny,1:nz]*dudx[2:nx,1:ny,1:nz] +                             0.5*(0.5*(v_o[1:nx-1,1:ny,1:nz] + v_o[1:nx-1,2:ny+1,1:nz])+0.5*(v_o[2:nx,1:ny,1:nz] + v_o[2:nx,2:ny+1,1:nz]))*dudy[2:nx,1:ny,1:nz] +                             0.5*(0.5*(w_o[1:nx-1,1:ny,1:nz] + w_o[1:nx-1,1:ny,2:nz+1])+0.5*(w_o[2:nx,1:ny,1:nz] + w_o[2:nx,1:ny,2:nz+1]))*dudz[2:nx,1:ny,1:nz] 
        
        H_vo[1:nx,2:ny,1:nz] = 0.5*(0.5*(u_o[1:nx,1:ny-1,1:nz]+u_o[2:nx+1,1:ny-1,1:nz]) + 0.5*(u_o[1:nx,2:ny,1:nz]+u_o[2:nx+1,2:ny,1:nz]))*dvdx[1:nx,2:ny,1:nz] +                             v_o[1:nx,2:ny,1:nz]*dvdy[1:nx,2:ny,1:nz] +                             0.5*(0.5*(w_o[1:nx,1:ny-1,1:nz] + w_o[1:nx,1:ny-1,2:nz+1])+0.5*(w_o[1:nx,2:ny,1:nz]+w_o[1:nx,2:ny,2:nz+1]))*dvdz[1:nx,2:ny,1:nz]
        
        H_wo[1:nx,1:ny,2:nz] = 0.5*(0.5*(u_o[1:nx,1:ny,1:nz-1]+u_o[2:nx+1,1:ny,1:nz-1]) + 0.5*(u_o[1:nx,1:ny,2:nz]+ u_o[2:nx+1,1:ny,2:nz]))*dwdx[1:nx,1:ny,2:nz] +                             0.5*(0.5*(v_o[1:nx,1:ny,1:nz-1] + v_o[1:nx,2:ny+1,1:nz-1]) + 0.5*(v_o[1:nx,1:ny,2:nz] + v_o[1:nx,2:ny+1,2:nz]))*dwdy[1:nx,1:ny,2:nz] +                             w_o[1:nx,1:ny,2:nz]*dwdz[1:nx,1:ny,2:nz]
                            
        
        dudx_2[2:nx,1:ny,1:nz] = ((u_o[3:nx+1,1:ny,1:nz] - u_o[2:nx,1:ny,1:nz])/dx - (u_o[2:nx,1:ny,1:nz] - u_o[1:nx-1,1:ny,1:nz])/dx)/dx
        dudy_2[1:nx,1:ny,1:nz] = ((u_o[1:nx,2:ny+1,1:nz] - u_o[1:nx,1:ny,1:nz])/dy - (u_o[1:nx,1:ny,1:nz] - u_o[1:nx,0:ny-1,1:nz])/dy)/dy
        dudz_2[1:nx,1:ny,1:nz] = ((u_o[1:nx,1:ny,2:nz+1] - u_o[1:nx,1:ny,1:nz])/dz - (u_o[1:nx,1:ny,1:nz] - u_o[1:nx,1:ny,0:nz-1])/dz)/dz
        
        dvdx_2[1:nx,1:ny,1:nz] = ((v_o[2:nx+1,1:ny,1:nz] - v_o[1:nx,1:ny,1:nz])/dx - (v_o[1:nx,1:ny,1:nz] - v_o[0:nx-1,1:ny,1:nz])/dx)/dx
        dvdy_2[1:nx,2:ny,1:nz] = ((v_o[1:nx,3:ny+1,1:nz] - v_o[1:nx,2:ny,1:nz])/dy - (v_o[1:nx,2:ny,1:nz] - v_o[1:nx,1:ny-1,1:nz])/dy)/dy
        dvdz_2[1:nx,1:ny,1:nz] = ((v_o[1:nx,1:ny,2:nz+1] - v_o[1:nx,1:ny,1:nz])/dy - (v_o[1:nx,1:ny,1:nz] - v_o[1:nx,1:ny,0:nz-1])/dz)/dz

        dwdx_2[1:nx,1:ny,1:nz] = ((w_o[2:nx+1,1:ny,1:nz] - w_o[1:nx,1:ny,1:nz])/dx - (w_o[1:nx,1:ny,1:nz] - w_o[0:nx-1,1:ny,1:nz])/dx)/dx
        dwdy_2[1:nx,1:ny,1:nz] = ((w_o[1:nx,2:ny+1,1:nz] - w_o[1:nx,1:ny,1:nz])/dy - (w_o[1:nx,1:ny,1:nz] - w_o[1:nx,0:ny-1,1:nz])/dy)/dy
        dwdz_2[1:nx,1:ny,2:nz] = ((w_o[1:nx,1:ny,3:nz+1] - w_o[1:nx,1:ny,2:nz])/dy - (w_o[1:nx,1:ny,2:nz] - w_o[1:nx,1:ny,1:nz-1])/dz)/dz
        #Debugging 여기에 0.5가 곱해져 있었음
        
        if t_real ==0: 
            H_up[2:nx,1:ny,1:nz] = H_uo[2:nx,1:ny,1:nz]; H_vp[1:nx,2:ny,1:nz] = H_vo[1:nx,2:ny,1:nz]
            H_wp[1:nx,1:ny,2:nz] = H_wo[1:nx,1:ny,2:nz] #only first step
            
        #d coefficient (Boundary condition은 따로 처리해줘야 한다)
        d_u[:] = dt/(Re)*(dudx_2[2:nx,1:ny,1:nz] + dudy_2[2:nx,1:ny,1:nz] + dudz_2[2:nx,1:ny,1:nz]) - dt/2*(3*H_uo[2:nx,1:ny,1:nz] - H_up[2:nx,1:ny,1:nz])
        d_v[:] = dt/(Re)*(dvdx_2[1:nx,2:ny,1:nz] + dvdy_2[1:nx,2:ny,1:nz] + dvdz_2[1:nx,2:ny,1:nz]) - dt/2*(3*H_vo[1:nx,2:ny,1:nz] - H_vp[1:nx,2:ny,1:nz])
        d_w[:] = dt/(Re)*(dwdx_2[1:nx,1:ny,2:nz] + dwdy_2[1:nx,1:ny,2:nz] + dwdz_2[1:nx,1:ny,2:nz]) - dt/2*(3*H_wo[1:nx,1:ny,2:nz] - H_wp[1:nx,1:ny,2:nz])
        
        #------------- u TDMA  ---------------#
        #x-direction TDMA for u velocity
        a_x[2:nx] = -mu; b_x[2:nx] = 1+2*mu; c_x[2:nx] = -mu
        
        x_d=1; y_d=0; z_d=0
        g[2:nx,1:ny,1:nz,:] = TDMA(a_x[2:nx],b_x[2:nx],c_x[2:nx],d_u,g[2:nx,1:ny,1:nz],factor_x[2:nx],nx-2,ny-1,nz-1,x_d,y_d,z_d)

        #y-direction TDMA for u velocity
        a_y[:,1:ny] = -mu; b_y[:,1:ny] = 1+2*mu; c_y[:,1:ny] = -mu
        b_y[:,1] = (1+3*mu); b_y[:,ny-1] = (1+3*mu)
        
        x_d = 0; y_d = 1; z_d = 0
        h[2:nx,1:ny,1:nz,:] = TDMA(a_y[:,1:ny],b_y[:,1:ny],c_y[:,1:ny],g[2:nx,1:ny,1:nz],h[2:nx,1:ny,1:nz],factor_y[:,1:ny],nx-2,ny-1,nz-1,x_d,y_d,z_d)
 
        #z-direction TDMA for u velocity
        a_z[:,:,1:nz] = -mu; b_z[:,:,1:nz] = 1+2*mu; c_z[:,:,1:nz] = -mu
        b_z[:,:,1] = (1+3*mu); b_z[:,:,nz-1] = (1+3*mu)
        
        x_d = 0; y_d = 0; z_d = 1
        u_mid[2:nx,1:ny,1:nz,:] = TDMA(a_z[:,:,1:nz],b_z[:,:,1:nz],c_z[:,:,1:nz],h[2:nx,1:ny,1:nz],q[2:nx,1:ny,1:nz],factor_z[:,:,1:nz],nx-2,ny-1,nz-1,x_d,y_d,z_d)
        
        #------------ v TDMA -----------------#
        #x-direction TDMA for v velocity
        a_x[1:nx] = -mu; b_x[1:nx] = 1+2*mu; c_x[1:nx] = -mu
        b_x[1,:] = (1+3*mu); b_x[nx-1,:] = (1+3*mu)
        
        x_d=1; y_d=0; z_d=0
        g[1:nx,2:ny,1:nz,:] = TDMA(a_x[1:nx],b_x[1:nx],c_x[1:nx],d_v,g[1:nx,2:ny,1:nz],factor_x[1:nx],nx-1,ny-2,nz-1,x_d,y_d,z_d)

        #y-direction TDMA for v velocity
        a_y[:,2:ny] = -mu; b_y[:,2:ny] = 1+2*mu; c_y[:,2:ny] = -mu
        
        x_d=0; y_d=1; z_d=0
        h[1:nx,2:ny,1:nz,:] = TDMA(a_y[:,2:ny],b_y[:,2:ny],c_y[:,2:ny],g[1:nx,2:ny,1:nz],h[1:nx,2:ny,1:nz],factor_y[:,2:ny],nx-1,ny-2,nz-1,x_d,y_d,z_d)

        #z-direction TDMA for v velocity
        a_z[:,:,1:nz] = -mu; b_z[:,:,1:nz] = 1+2*mu; c_z[:,:,1:nz] = -mu
        b_z[:,:,1] = (1+3*mu); b_z[:,:,nz-1] = (1+3*mu)
        
        x_d=0; y_d=0; z_d=1
        
        v_mid[1:nx,2:ny,1:nz,:] = TDMA(a_z[:,:,1:nz],b_z[:,:,1:nz],c_z[:,:,1:nz],h[1:nx,2:ny,1:nz],q[1:nx,2:ny,1:nz],factor_z[:,:,1:nz],nx-1,ny-2,nz-1,x_d,y_d,z_d)
        
        #------------ w TDMA -----------------#
        #x-direction TDMA for w velocity
        a_x[1:nx] = -mu; b_x[1:nx] = 1+2*mu; c_x[1:nx] = -mu
        b_x[1,:] = (1+3*mu); b_x[nx-1,:] = (1+3*mu)
        
        x_d=1; y_d=0; z_d=0
        g[1:nx,1:ny,2:nz,:] = TDMA(a_x[1:nx],b_x[1:nx],c_x[1:nx],d_w,g[1:nx,1:ny,2:nz],factor_x[1:nx],nx-1,ny-1,nz-2,x_d,y_d,z_d)

        #y-direction TDMA for w velocity
        a_y[:,1:ny] = -mu; b_y[:,1:ny] = 1+2*mu; c_y[:,1:ny] = -mu
        b_y[:,1] = (1+3*mu); b_y[:,ny-1] = (1+3*mu)        
        
        x_d=0; y_d=1; z_d=0
        h[1:nx,1:ny,2:nz,:] = TDMA(a_y[:,1:ny],b_y[:,1:ny],c_y[:,1:ny],g[1:nx,1:ny,2:nz],h[1:nx,1:ny,2:nz],factor_y[:,1:ny],nx-1,ny-1,nz-2,x_d,y_d,z_d)

        #z-direction TDMA for w velocity
        a_z[:,:,2:nz] = -mu; b_z[:,:,2:nz] = 1+2*mu; c_z[:,:,2:nz] = -mu

        x_d=0; y_d=0; z_d=1
        w_mid[1:nx,1:ny,2:nz,:] = TDMA(a_z[:,:,2:nz],b_z[:,:,2:nz],c_z[:,:,2:nz],h[1:nx,1:ny,2:nz],q[1:nx,1:ny,2:nz],factor_z[:,:,2:nz],nx-1,ny-1,nz-2,x_d,y_d,z_d)
                
        #update intermeidate velocity
        u_mid[2:nx,1:ny,1:nz] = u_o[2:nx,1:ny,1:nz] + u_mid[2:nx,1:ny,1:nz]
        v_mid[1:nx,2:ny,1:nz] = v_o[1:nx,2:ny,1:nz] + v_mid[1:nx,2:ny,1:nz]
        w_mid[1:nx,1:ny,2:nz] = w_o[1:nx,1:ny,2:nz] + w_mid[1:nx,1:ny,2:nz]
        
        #-------------- Poisson equation ----------------#
        p_o, steps = Poisson(p_n,p_o,u_mid,v_mid,w_mid,dt,nx,ny,nz,dx,dy,dz,epsil_p)
        #Poisson equation에서 update된 p_o를 받아와야 하는데 p_n을 받아와서 터진적이 있음
        
        #Correction pressure term
        u_n[2:nx,1:ny,1:nz,:] = u_mid[2:nx,1:ny,1:nz,:] - dt*(p_o[2:nx,1:ny,1:nz] - p_o[1:nx-1,1:ny,1:nz])/dx #u_update
        v_n[1:nx,2:ny,1:nz,:] = v_mid[1:nx,2:ny,1:nz,:] - dt*(p_o[1:nx,2:ny,1:nz] - p_o[1:nx,1:ny-1,1:nz])/dy #v_update
        w_n[1:nx,1:ny,2:nz,:] = w_mid[1:nx,1:ny,2:nz,:] - dt*(p_o[1:nx,1:ny,2:nz] - p_o[1:nx,1:ny,1:nz-1])/dz #w_update
        
        print("Residual")
        print(np.sum(abs(u_n[1:nx,1:ny,1:nz] - u_o[1:nx,1:ny,1:nz]) + abs(v_n[1:nx,1:ny,1:nz] - v_o[1:nx,1:ny,1:nz])                      + abs(w_n[1:nx,1:ny,1:nz] - w_o[1:nx,1:ny,1:nz]))/((nx-1)*(ny-1)*(nz-1)))
        
        Residual = (np.sum(abs(u_n[1:nx,1:ny,1:nz] - u_o[1:nx,1:ny,1:nz]) + abs(v_n[1:nx,1:ny,1:nz] - v_o[1:nx,1:ny,1:nz])                      + abs(w_n[1:nx,1:ny,1:nz] - w_o[1:nx,1:ny,1:nz]))/((nx-1)*(ny-1)*(nz-1)))
        
        if (Residual < epsil_u):
            break
        
        u_o[1:nx,1:ny,1:nz] = u_n[1:nx,1:ny,1:nz]; v_o[1:nx,1:ny,1:nz] = v_n[1:nx,1:ny,1:nz]
        w_o[1:nx,1:ny,1:nz] = w_n[1:nx,1:ny,1:nz]
        
        
        
        if t_real %50 ==0:
            
            #Staggered to Collocated grids
            u_c[0:nx,0:ny,0:nz,0:1] = 0.5*(0.5*(u_n[1:nx+1,0:ny,0:nz] + u_n[1:nx+1,1:ny+1,0:nz])+0.5*(u_n[1:nx+1,0:ny,1:nz+1] + u_n[1:nx+1,1:ny+1,1:nz+1]))
                
            u_c[0:nx,0:ny,0:nz,1:2] = 0.5*(0.5*(v_n[0:nx,1:ny+1,0:nz] + v_n[1:nx+1,1:ny+1,0:nz])+0.5*(v_n[0:nx,1:ny+1,1:nz+1] + v_n[1:nx+1,1:ny+1,1:nz+1]))
            
            u_c[0:nx,0:ny,0:nz,2:3] = 0.5*(0.5*(w_n[0:nx,0:ny,1:nz+1] + w_n[1:nx+1,0:ny,1:nz+1])+0.5*(w_n[0:nx,1:ny+1,1:nz+1] + w_n[1:nx+1,1:ny+1,1:nz+1]))
            
            #get vorticity here #z방향 위에 + 아래 평균
            vorticity[:,:,:,0:1] = 0.5*(0.5*((u_c[1:nx,1:ny,0:nz-1,2:3]-u_c[1:nx,0:ny-1,0:nz-1,2:3])/dy + (u_c[0:nx-1,1:ny,0:nz-1,2:3]-u_c[0:nx-1,0:ny-1,0:nz-1,2:3])/dy)+                  0.5*((u_c[1:nx,1:ny,1:nz,2:3]-u_c[1:nx,0:ny-1,1:nz,2:3])/dy + (u_c[0:nx-1,1:ny,1:nz,2:3]-u_c[0:nx-1,0:ny-1,1:nz,2:3])/dy))-            0.5*(0.5*((u_c[0:nx-1,1:ny,1:nz,1:2] - u_c[0:nx-1,1:ny,0:nz-1,1:2])/dz + (u_c[0:nx-1,0:ny-1,1:nz,1:2] - u_c[0:nx-1,0:ny-1,0:nz-1,1:2])/dz)+                 0.5*((u_c[1:nx,1:ny,1:nz,1:2] - u_c[1:nx,1:ny,0:nz-1,1:2])/dz + (u_c[1:nx,0:ny-1,1:nz,1:2] - u_c[1:nx,0:ny-1,0:nz-1,1:2])/dz))

            
            vorticity[:,:,:,1:2] = 0.5*(0.5*((u_c[0:nx-1,1:ny,1:nz,0:1] - u_c[0:nx-1,1:ny,0:nz-1,0:1])/dz + (u_c[0:nx-1,0:ny-1,1:nz,0:1] - u_c[0:nx-1,0:ny-1,0:nz-1,0:1])/dz)+                 0.5*((u_c[1:nx,1:ny,1:nz,0:1] - u_c[1:nx,1:ny,0:nz-1,0:1])/dz + (u_c[1:nx,0:ny-1,1:nz,0:1] - u_c[1:nx,0:ny-1,0:nz-1,0:1])/dz))-            0.5*(0.5*((u_c[1:nx,1:ny,0:nz-1,2:3] - u_c[0:nx-1,1:ny,0:nz-1,2:3])/dx + (u_c[1:nx,0:ny-1,0:nz-1,2:3] - u_c[0:nx-1,0:ny-1,0:nz-1,2:3])/dx)+                 0.5*((u_c[1:nx,1:ny,1:nz,2:3] - u_c[0:nx-1,1:ny,1:nz,2:3])/dx + (u_c[1:nx,0:ny-1,1:nz,2:3] - u_c[0:nx-1,0:ny-1,1:nz,2:3])/dx))
                
                
            vorticity[:,:,:,2:3] = 0.5*(0.5*((u_c[1:nx,1:ny,0:nz-1,1:2] - u_c[0:nx-1,1:ny,0:nz-1,1:2])/dx + (u_c[1:nx,0:ny-1,0:nz-1,1:2] - u_c[0:nx-1,0:ny-1,0:nz-1,1:2])/dx)+                 0.5*((u_c[1:nx,1:ny,1:nz,1:2] - u_c[0:nx-1,1:ny,1:nz,1:2])/dx + (u_c[1:nx,0:ny-1,1:nz,1:2] - u_c[0:nx-1,0:ny-1,1:nz,1:2])/dx))-            0.5*(0.5*((u_c[1:nx,1:ny,0:nz-1,0:1]-u_c[1:nx,0:ny-1,0:nz-1,0:1])/dy + (u_c[0:nx-1,1:ny,0:nz-1,0:1]-u_c[0:nx-1,0:ny-1,0:nz-1,0:1])/dy)+                  0.5*((u_c[1:nx,1:ny,1:nz,0:1]-u_c[1:nx,0:ny-1,1:nz,0:1])/dy + (u_c[0:nx-1,1:ny,1:nz,0:1]-u_c[0:nx-1,0:ny-1,1:nz,0:1])/dy))

            
            #print contour stream line & velocity field
            video_3D = open('3D u_fina_Video.plt', 'a')

            if t_real==0: video_3D.write('VARIABLES = "x","y","z","u","v","w"\n')
            video_3D.write('Zone T="HIT%d"\n'%0)
            video_3D.write('I=%d J=%d K=%d\n' %(nx-1, ny-1, nz-1))
    
            for k in range(1, nz):
                for j in range(1, ny):
                    for i in range(1, nx):
                        video_3D.write('%.5f %.5f %.5f %.5f %.5f %.5f\n' %(i,j,k, u_o[i,j,k,0], v_o[i,j,k,0], w_o[i,j,k,0]))

    
            video_3D.close()    
            
            #vorticity field
            Vorticity_video = open('3D vorticity.plt', 'a')

            if t_real==0: Vorticity_video.write('VARIABLES = "x","y","z","vor_x","vor_y","vor_z"\n')
            Vorticity_video.write('Zone T="HIT%d"\n'%0)
            Vorticity_video.write('I=%d J=%d K=%d\n' %(nx-1, ny-1, nz-1))
    
            for k in range(nz-1):
                for j in range(ny-1):
                    for i in range(nx-1):
                        Vorticity_video.write('%.5f %.5f %.5f %.5f %.5f %.5f\n' %(i,j,k, vorticity[i,j,k,0], vorticity[i,j,k,1], vorticity[i,j,k,2]))

    
            Vorticity_video.close()            

        
        t_real += 1
    
    duration = timer() - start 
        
    #Boundary condition
    u_o[1:nx,ny,1:nz,:] = (2 - u_o[1:nx,ny-1,1:nz,:]); u_o[1:nx,0,1:nz,:] = - u_o[1:nx,1,1:nz,:]
    u_o[1:nx,1:ny,nz,:] = -u_o[1:nx,1:ny,nz-1,:]; u_o[1:nx,1:ny,0,:] = - u_o[1:nx,1:ny,1,:]
        
    v_o[nx,1:ny,1:nz,:] = - v_o[nx-1,1:ny,1:nz,:]; v_o[0,1:ny,1:nz,:] = - v_o[1,1:ny,1:nz,:]
    v_o[1:nx,1:ny,nz,:] = - v_o[1:nx,1:ny,nz-1,:]; v_o[1:nx,1:ny,0,:] = - v_o[1:nx,1:ny,1,:]

    w_o[nx,1:ny,1:nz,:] = - w_o[nx-1,1:ny,1:nz,:]; w_o[0,1:ny,1:nz,:] = - w_o[1,1:ny,1:nz,:]
    w_o[1:nx,ny,1:nz,:] = - w_o[1:nx,ny-1,1:nz,:]; w_o[1:nx,0,1:nz,:] = - w_o[1:nx,1,1:nz,:]
    
    #Staggered grids to Collocated grids
    u_c[0:nx,0:ny,0:nz,0:1] = (u_n[1:nx+1,0:ny,0:nz] + u_n[1:nx+1,1:ny+1,1:nz+1])/2
    u_c[0:nx,0:ny,0:nz,1:2] = (v_n[0:nx,1:ny+1,0:nz] + v_n[1:nx+1,1:ny+1,1:nz+1])/2  
    u_c[0:nx,0:ny,0:nz,2:3] = (w_n[0:nx,0:nz,1:ny+1] + w_n[1:nx+1,1:ny+1,1:nz+1])/2
            
    
    #print contour stream line & velocity field
    contour = open('3D u_fina.plt', 'w')

    contour.write('VARIABLES = "x","y","z","u","v","w"\n')
    contour.write('Zone T="HIT%d"\n'%0)
    contour.write('I=%d J=%d K=%d\n' %(nx-1, ny-1, nz-1))
    
    for k in range(1, nz):
        for j in range(1, ny):
            for i in range(1, nx):
                contour.write('%.5f %.5f %.5f %.5f %.5f %.5f\n' %(i,j,k, u_o[i,j,k,0], v_o[i,j,k,0], w_o[i,j,k,0]))
    
    contour.close()
    
    print("end")
    print(duration)
    print(t_real)
            
        
if __name__ == "__main__":
        
    main()
        
    print("All process is finished!")
    

