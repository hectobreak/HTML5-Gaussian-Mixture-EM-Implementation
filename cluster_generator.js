
function erf(x, k){
  let c   = x / k / 6;
  let C4  = Math.exp(- 0.25 * x * x / k / k);
  let C   = C4 * C4 * C4 * C4;
  let Cn  = 1;
  let Cn2 = 1;
  let sum = 0;
  for(let n = 1; n < k; ++n){
    Cn  *= 1 / C;
    Cn2 *= 1 / C / Cn / Cn;
    sum += c * (Cn2 * Cn * Cn * C + 4 * Cn2 * Cn * C4 + Cn2);
  }  
  return sum * 2 / Math.sqrt(Math.PI);
}

function gaussianSample(){
	let randomVal = Math.random();
	let gInt = (x) => 0.5 + erf(x, 200)/2;
	let ini = -5, end = 5;
	for(let i = 0; i < 50; ++i){
		let mid = (ini + end)/2;
		if(gInt(mid) > randomVal) end = mid;
		else ini = mid;
	}
	return (ini + end) / 2;
}

function rotatePoint(point, angle){
	let x = point[0], y = point[1];
	point[0] =   x * Math.cos(angle) + y * Math.sin(angle);
	point[1] = - x * Math.sin(angle) + y * Math.cos(angle);
}

function generateClusters(nclusters, clusterSizeMin, clusterSizeMax, sizeMin, sizeMax) {
	let points = new Array;
	
	let canvasWidth = 1000, canvasHeight = 600;
	// 1000 x 600
	for(let i = 0; i < nclusters; ++i){
		let npoints = Math.floor((clusterSizeMax - clusterSizeMin) * Math.random() + clusterSizeMin);
		let sx =  Math.floor((sizeMax - sizeMin) * Math.random() + sizeMin);
		let sy =  Math.floor((sizeMax - sizeMin) * Math.random() + sizeMin);
		let rx = Math.random() * 2 * Math.PI;
		let tx = Math.floor(Math.random() * 800) - 400;
		let ty = Math.floor(Math.random() * 400) - 200;
		for(let n = 0; n < npoints; ++n){
			  let newPoint = [ gaussianSample() * sx, gaussianSample() * sy ];
			  rotatePoint(newPoint, rx);
			  newPoint[0] += tx;
			  newPoint[1] += ty;
			  points.push(new Vector(newPoint));
		}
		console.log(tx, ty);
	}
	return points;
}
