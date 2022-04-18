class Vector {
	constructor( arr ){
		if(arr instanceof Vector){
			this.size = arr.size;
			this.vec = arr.vec.slice();
		} else {
			this.size = arr.length;
			this.vec = arr;
		}
	}
	
	dot( vec ){
		if( !(vec instanceof Vector) ) throw new Error("The input must be a vector");
		if( vec.size !== this.size ) throw new Error("The sizes have to be the same");
		let sum = 0; 
		for(let i = 0; i < this.size; ++i) sum += this.vec[i] * vec.vec[i];
		return sum;
	}
	
	add( vec ){
		if( !(vec instanceof Vector) ) throw new Error("The input must be a vector");
		if( vec.size !== this.size ) throw new Error("The sizes have to be the same");
		let sum = new Array(vec.size); 
		for(let i = 0; i < this.size; ++i) sum[i] = this.vec[i] + vec.vec[i];
		return new Vector(sum);
	}
	
	sub( vec ){
		if( !(vec instanceof Vector) ) throw new Error("The input must be a vector");
		if( vec.size !== this.size ) throw new Error("The sizes have to be the same");
		let sum = new Array(vec.size); 
		for(let i = 0; i < this.size; ++i) sum[i] = this.vec[i] - vec.vec[i];
		return new Vector(sum);
	}
	
	constructMatrix( vec ) {
		if( !(vec instanceof Vector) ) throw new Error("The input must be a vector");
		if( vec.size !== this.size ) throw new Error("The sizes have to be the same");
		let m = new Array(vec.size * vec.size);
		for(let i = 0; i < vec.size; ++i){
			for(let j = 0; j < vec.size; ++j)
				m[i + j * vec.size] = this.vec[i] * vec.vec[j];
		}
		return new Matrix(m);
	}
	
	mult( val ){
		if( typeof val !== "number" ) {
			throw new Error("The input must be a number");
		}
		let v = new Array( this.size );
		for(let i = 0; i < this.size; ++i) v[i] = this.vec[i] * val;
		return new Vector(v);
	}
}

class Matrix {
	constructor( arr ){
		if(typeof arr === "number") {
			this.size = arr;
			this.matrix = new Array( arr * arr );
			for(let i = 0; i < arr * arr; ++i) this.matrix[i] = 0;
		} else if( arr instanceof Matrix ) {
			this.size = arr.size;
			this.matrix = arr.matrix.slice();
			if(arr.triangular !== undefined) 
				this.triangular = arr.triangular.slice();
			if(arr.triangularRes !== undefined) 
				this.triangularRes = arr.triangularRes.slice();
			if(arr.numberFlips !== undefined) 
				this.numberFlips = arr.numberFlips;
			if(arr.det !== undefined) 
				this.det = arr.det;
			if(arr.inv !== undefined) 
				this.inv = arr.inv.slice();
		} else {
			let n0, n, m = 2;
			do {
				[n0, n, m] = [n, m, Math.floor( 0.5 * (m + arr.length / m) )];
			} while(n !== m && n0 !== m);
			if(n0 !== undefined) n = Math.min(n0, n, m);
			while(n*n < arr.length) n++;
			if(n*n !== arr.length) {
				throw new Error("The matrix is not a square matrix");
			}
			this.size = n;
			this.matrix = arr;
		}
	}
	
	makeTriangular(){
		if(this.triangular !== undefined) return new Matrix(this.triangular);
		
		this.triangular    = new Array(this.matrix.length);
		this.triangularRes = new Array(this.matrix.length);
		this.numberFlips   = 0;
		
		for(let i = 0; i < this.matrix.length; ++i) {
			this.triangular[i] = this.matrix[i];
			this.triangularRes[i] = 0;
		}
		for(let i = 0; i < this.size; ++i)
			this.triangularRes[i + this.size*i] = 1;
		
		for(let i = 0; i < this.size; ++i){
			// Ensure pivot is non-zero
			let p = i;
			while( p < this.size && this.triangular[i + p*this.size] === 0 ) ++p;
			if(p === this.size) continue;
			if(p !== i) {
				this.numberFlips++;
				for(let j = i; j < this.size; ++j){
					[ this.triangular[j+i*this.size], this.triangular[j+p*this.size] ]
					= [ this.triangular[j+p*this.size], this.triangular[j+i*this.size] ];
				}
				for(let j = 0; j < this.size; ++j){
					[ this.triangularRes[j+i*this.size], this.triangularRes[j+p*this.size] ]
					= [ this.triangularRes[j+p*this.size], this.triangularRes[j+i*this.size] ];
				}
			}
			
			// Clear column i
			for(let p = i + 1; p < this.size; ++p){
				let factor = this.triangular[i+p*this.size] / this.triangular[i+i*this.size];
				for(let j = i; j < this.size; ++j){
					this.triangular[j+p*this.size] -= factor * this.triangular[j+i*this.size];
				}
				for(let j = 0; j < this.size; ++j){
					this.triangularRes[j+p*this.size] -= factor * this.triangularRes[j+i*this.size];
				}
			}
		}
		return new Matrix(this.triangular);
	}
	
	determinant(){
		if(this.det !== undefined) return this.det;
		if(this.triangular === undefined) this.makeTriangular();
		this.det = 1 - 2 * (this.numberFlips % 2);
		for(let i = 0; i < this.size; ++i) this.det *= this.triangular[i + this.size*i];
		return this.det;
	}
	
	inverse(){
		if(this.inv !== undefined) {
			let m1 = new Matrix(this.inv);
			m1.inv = this.matrix;
			return m1;
		}
		if(this.det === undefined) this.determinant();
		if(this.det === 0) throw new Error("The matrix is not invertible");
		
		this.inv = this.triangularRes.slice();
		let tmp = this.triangular.slice();
		
		for(let i = this.size - 1; i >= 0; --i){
			// Clear column i
			for(let p = i - 1; p >= 0; --p){
				let factor = tmp[i+p*this.size] / tmp[i+i*this.size];
				for(let j = i; j < this.size; ++j){
					tmp[j+p*this.size] -= factor * tmp[j+i*this.size];
				}
				for(let j = 0; j < this.size; ++j){
					this.inv[j+p*this.size] -= factor * this.inv[j+i*this.size];
				}
			}
		}
		
		for(let i = 0; i < this.size; ++i){
			let div = tmp[i + i*this.size];
			for(let j = 0; j < this.size; ++j)
				this.inv[j+i*this.size] /= div;
		}
		let m1 = new Matrix(this.inv);
		m1.inv = this.matrix;
		return m1;
	}
	
	transpose(){
		if(this.trans !== undefined) {
			let m1 = new Matrix(this.trans);
			m1.trans = this.matrix;
			return m1;
		}
		this.trans = new Array(this.matrix.length);
		for(let i = 0; i < this.size; ++i){
			for(let j = 0; j < this.size; ++j){
				this.trans[i+j*this.size] = this.matrix[j+i*this.size];
			}
		}
		
		let m1 = new Matrix(this.trans);
		m1.trans = this.matrix;
		return m1;
	}
	
	equals( matr ){
		if( !( val instanceof Matrix )) return false;
		return this.matrix.join(";") === matr.matrix.join(";");
	}
	
	isSymmetric(){
		return this.equals( this.transpose() );
	}
	
	isPositiveDefinite(){
		if(!this.isSymmetric()) return false;
		/* Find out how to determine this */
	}
	
	mult( val ){
		if( val instanceof Matrix ){
			if(val.size !== this.size) throw new Error("Matrices must be of the same size");
			let arr = new Array(this.size * this.size);
			for(let i = 0; i < arr.length; ++i) arr[i] = 0;
			
			for(let i = 0; i < val.size; ++i){
				for(let j = 0; j < val.size; ++j){
					for(let k = 0; k < val.size; ++k){
						arr[i+j*val.size] += this.matrix[k+ j*val.size] * val.matrix[i + k*val.size];
					}
				}
			}
			return new Matrix(arr);
		} else if( val instanceof Vector ){
			if(val.size !== this.size) throw new Error("Vectors must be of the same size");
			let arr = new Array(this.size);
			for(let i = 0; i < arr.length; ++i) arr[i] = 0;
			
			for(let i = 0; i < val.size; ++i){
				for(let j = 0; j < val.size; ++j){
					arr[i] += this.matrix[j+ i*val.size] * val.vec[j];
				}
			}
			return new Vector(arr);
		} else if( typeof val === 'number' ) {
			let m = new Array(this.matrix.length);
			for(let i = 0; i < this.matrix.length; ++i)
				m[i] = this.matrix[i] * val;
			return new Matrix(m);
		} else throw new Error("Invalid input");
	}
	
	add( val ){
		if( val instanceof Matrix ){
			if(val.size !== this.size) throw new Error("Matrices must be of the same size");
			let arr = new Array(this.size * this.size);
			
			for(let i = 0; i < val.matrix.length; ++i){
				arr[i] = this.matrix[i] + val.matrix[i];
			}
			return new Matrix(arr);
		} else throw new Error("Invalid input");
	}
}
