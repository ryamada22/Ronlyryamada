
tmp.out <- my.mesh.back.euler(vertices,faces.v,step=0.05)
my.mesh.tri.plot(tmp.out,faces.v)
tmp.out2 <- my.mesh.back.euler(tmp.out,faces.v,step=0.05)
my.mesh.tri.plot(tmp.out2,faces.v)



my.mesh.back.euler <- function(vertices,faces.v,step = 0.01,eps = 10^-5,max.iter = 50){
	# 曲率,rho と法線ベクトル,N
	rho.out <- my.curvature.cot(vertices,faces.v)
	rho.N <- rho.out[[1]]
	new.vertices <- vertices - 2 * rho.N * step
	for(i in 1:max.iter){
		# new.verticesでのrho=new.rho,N=new.N
		rho.out <- my.curvature.cot(new.vertices,faces.v)
		new.rho.N <- rho.out[[1]]
		tmp.vertices <- vertices - 2 * new.rho.N * step
		diff.vertices <- tmp.vertices - new.vertices
		if(max(Norm(diff.vertices) < eps)){
			break
		}
		new.vertices = tmp.vertices		
	}
	new.vertices
}

my.mesh.tri.plot(vertices,faces.v)
tmp.out <- my.mesh.back.euler(vertices,faces.v,step=0.005,eps=10^(-10),max.iter=200)
for(i in 1:20){
	tmp.out <- my.mesh.back.euler(tmp.out,faces.v,step=0.005)
	my.mesh.tri.plot(tmp.out,faces.v)
}
my.mesh.tri.plot(tmp.out,faces.v)
