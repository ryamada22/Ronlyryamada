my.mesh.back.euler <- function(vertices,faces.v,step = 0.01,eps = 10^-5,max.iter = 50){
	# 曲率,rho と法線ベクトル,N
	rho.out <- my.conformal.rho(vertices,faces.v)
	rho.N <- rho.out[[3]] * Im(rho.out[[2]])
	new.vertices <- vertices + 2 * rho.N * step
	for(i in 1:max.iter){
		# new.verticesでのrho=new.rho,N=new.N
		rho.out <- my.conformal.rho(new.vertices,faces.v)
		new.rho.N <- rho.out[[3]] * Im(rho.out[[2]])
		tmp.vertices <- vertices + 2 * new.rho.N * step
		diff.vertices <- tmp.vertices - new.vertices
		if(max(Norm(diff.vertices) < epx)){
			break
		}
		new.vertices = tmp.vertices		
	}
	new.vertices
}