class Gaussian {
	constructor( dim, sigma ){
		if(sigma !== undefined ){ 
			this.c = new Vector(dim);
			this.m = new Matrix(sigma);
			this.d = this.m.determinant();
			this.d2 = Math.pow( this.d, -0.5 );
			this.C = Math.pow( 2 * Math.PI, -dim.size / 2 );
			this.m1 = this.m.inverse();
		} else {
			let center = new Array(dim);
			let matrix = new Array(dim * dim);
			
			for(let i = 0; i < dim; ++i) {
				center[i] = (Math.random() - 0.5) * 600;
				for(let j = 0; j < dim; ++j)
					matrix[i + j*dim] = (i === j ? 10000 : 0);
			}
			this.c = new Vector(center);
			this.m = new Matrix(matrix);
			this.d = this.m.determinant();
			this.d2 = Math.pow( this.d, -0.5 );
			this.C = Math.pow( 2 * Math.PI, -dim / 2 );
			this.m1 = this.m.inverse();
		}
	}
	
	evaluate( x ){
		let dif = x.sub(this.c);
		return this.C * this.d2 * Math.exp( - 0.5 *  dif.dot( this.m1.mult( dif ) ) );
	}
	
}

function initialize( dataset, n ) {
	let Distributions = new Array(n);
	for(let i = 0; i < n; ++i){
		Distributions[i] = new Gaussian( dataset[0].size );
	}
	
	let Pi = new Array(n);
	for(let i = 0; i < dataset.length; ++i)
		Pi[i] = 1 / n;
	
	return [ Pi , Distributions ];
}

function ComputeGamma( dataset, Pi, Distributions ){
	let Gamma = new Array( dataset.length );
	for(let i = 0; i < dataset.length; ++i) {
		Gamma[i] = new Array( Distributions.length );
		let sum = 0;
		for(let c = 0; c < Distributions.length; ++c) {
			Gamma[i][c] = Distributions[c].evaluate( dataset[i] ) * Pi[c];
			sum += Gamma[i][c];
		}
		for(let c = 0; c < Distributions.length; ++c){
			Gamma[i][c] /= sum;
		}
	}
	return Gamma;
}

function recompute( dataset, Pi, Distributions, Gamma ){
	let mu = new Array( Distributions.length );
	let sigma = new Array( Distributions.length );
	let pi = new Array( Distributions.length );
	
	let gammasum = new Array( Distributions.length );
	
	for(let c = 0; c < Distributions.length; ++c){
		gammasum[c] = 0;
		for(let i = 0; i < dataset.length; ++i)
			gammasum[c] += Gamma[i][c];
	}
	
	for(let c = 0; c < Distributions.length; ++c){
		mu[c] = new Vector( Distributions[c].c.vec.map( x => 0 ) );
		sigma[c] = new Matrix(Distributions[c].m.size);
		pi[c] = gammasum[c] / Distributions.length;
		for(let i = 0; i < dataset.length; ++i){
			let xmu = dataset[i].sub( Distributions[c].c );
			mu[c] = mu[c].add( dataset[i].mult( Gamma[i][c] / gammasum[c] ) );
			sigma[c] = sigma[c].add( xmu.constructMatrix(xmu).mult( Gamma[i][c] / gammasum[c] ) );
		}
	}
	
	let distributions = new Array( Distributions.length );
	for(let c = 0; c < Distributions.length; ++c){
		distributions[c] = new Gaussian( mu[c], sigma[c] );
	}
	
	return [ pi, distributions ];
}

function computeL( dataset, Pi, Distributions ){
	let l = 0;
	for(let i = 0; i < dataset.length; ++i){
		let sum = 0;
		for(let c = 0; c < Distributions.length; ++c){
			sum += Pi[c] * Distributions[c].evaluate(dataset[i]);
		}
		l += Math.log( sum );
	}
	return l;
}

let td = 0.0001;

let allDistr = new Array;
let gammaDistr = new Array;

function ExpectationMaximization( dataset, n=2 ){
	let [ Pi, Distributions ] = initialize( dataset, n );
	let l0, l = Infinity;
	let idempotency = false;
	allDistr = new Array;
	gammaDistr = new Array;
	do {
		allDistr.push(Distributions);
		// Expectation
		let Gamma = ComputeGamma( dataset, Pi, Distributions );
		gammaDistr.push(Gamma);
		// Maximization
		[ Pi, Distributions ] = recompute( dataset, Pi, Distributions, Gamma );
		// Reevaluation
		l0 = l;
		l = computeL( dataset, Pi, Distributions );
		//console.log(l, l0, Math.abs( (l - l0) / l ));
	} while( Math.abs( (l - l0) ) >= td );
	allDistr.push(Distributions);
	gammaDistr.push(ComputeGamma( dataset, Pi, Distributions ));
	return Distributions;
}
