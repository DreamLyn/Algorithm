package com.locaris.algorithm.filter;

import com.locaris.algorithm.beans.Position;

import com.locaris.algorithm.residual.Residual;

/**
 * 定位滤波算法
 * @author Lyn
 */
public class LocationFilter {

	/**
	 * 残差值滤波，利用计算出的坐标值计算出据差，然后与真实据差进行比较。
	 * @param rdif
	 * @param devicePosition--->第一列和第二列分别为基站的XY坐标
	 * @param position
	 * @param residualValue
	 * @return
	 */
	public boolean residualFilter(double[][] rdif,double[][] devicePosition,Position position,double residualValue){
		double resudualValue= Residual.calculateResidual(rdif,devicePosition,position);
		if (resudualValue >= residualValue) {
			return false;
		}else{
			return true;
		}
	}

	/**
	 * 分别计算tempX和tempY的中值后组成Position
	 * @param tempX
	 * @param tempY
	 * @return
	 */
	public Position calculateCenter(double[] tempX,double[] tempY){
		int length=tempX.length;
		int centerLength= (int) Math.ceil(length/2.0);
		for (int i = 0; i < centerLength; i++) {
			for (int j = i + 1; j < length; j++) {
				if (tempX[i] > tempX[j]) {
					double tempDouble = tempX[i];
					tempX[i] = tempX[j];
					tempX[j] = tempDouble;
				}
				if (tempY[i] > tempY[j]) {
					double tempDouble = tempY[i];
					tempY[i] = tempY[j];
					tempY[j] = tempDouble;
				}
			}
		}
		Position position=new Position();
		position.setX(tempX[centerLength-1]);
		position.setY(tempY[centerLength-1]);
		return position;
	}
	/**
	 * 中值滤波
	 * @param positon
	 * @return
	 */
	public Position centerFilter(Position[] positon,int length){
		double tempX[] = new double[length];
		double tempY[] = new double[length];
		
		for (int count = 0; count < length; count++) {
			tempX[count] = positon[count].getX();
			tempY[count] = positon[count].getY();
		}
		return calculateCenter(tempX,tempY);

	}

	/**
	 * 获取绝对值计算后的中值
	 * @param positon
	 * @param uparam
	 * @return
	 */
	private Position centerFilterAbs(Position[] positon,Position uparam,int length){
		double tempX[] = new double[length];
		double tempY[] = new double[length];
		
		for (int count = 0; count < length; count++) {
			tempX[count] = Math.abs(positon[count].getX()-uparam.getX());
			tempY[count] = Math.abs(positon[count].getY()-uparam.getY());
		}
		return calculateCenter(tempX,tempY);
	}
	/**
	 * 高斯滤波
	 * @param param
	 * @param positions
	 */
	public Position gaussFilter(Position[] positions,int param,int length){
		Position position=new Position();
		double sumX=0,sumY=0,sumwX=0,sumwY=0;
		for(int count=0;count<length;count++){
			double tempX=1/(Math.sqrt(2*Math.PI)*param)*Math.exp(0-Math.pow(positions[(length-1)/2].getX()-positions[count].getX(),2)/2/(param*param));
			double tempY=1/(Math.sqrt(2*Math.PI)*param)*Math.exp(0-Math.pow(positions[(length-1)/2].getY()-positions[count].getY(),2)/2/(param*param));
			sumwX+=tempX;
			sumwY+=tempY;
			sumX+=positions[count].getX()*tempX;
			sumY+=positions[count].getY()*tempY;
		}
		double gaX=sumX/sumwX;
		double gaY=sumY/sumwY;
		Position positionM=centerFilter(positions,length);
		double uX=positionM.getX();
		double uY=positionM.getY();
		Position positionSigma=centerFilterAbs(positions, positionM,length);
		double sigmaX=1.483*positionSigma.getX();
		double sigmaY=1.483*positionSigma.getY();
		if((gaX<(uX-3*sigmaX))||(gaX>(uX+3*sigmaX))){
			gaX=uX;
		}
		
		if((gaY<(uY-3*sigmaY))||(gaY>(uY+3*sigmaY))){
			gaY=uY;
		}
		position.setX(gaX);
		position.setY(gaY);
		return position;
	}
	
	public Position sincFilter(Position[] positions,double bandwidth,double interval){
		Position position=new Position();
	
		return position;
	}
	
	
//	/**
//	 * 非参数滤波
//	 * @param positions
//	 * @return
//	 */
//	public Position nonparametricFilter(Position[] positions,int length,int index,int delay){
//		Position position=new Position();
//		double hParam=10;
//		int[] time0=new int[length];
//		int[] time1=new int[length];
//		double[][] x=new double[length][1];
//		double[][] y=new double[length][1];
////		double[] z=new double[length];
//		double[][] a=new double[length][2];
//		double[] k=new double[length];
//		double[][] w=new double[length][length];
//		for(int i=0;i<length;i++){
//			time0[i]=positions[i].getSequence()&0xff;
//		}
//		///////////
//		int ondex=(index==0)?(length-1):(index-1);
//		for(int i=0;i<length;i++){
//			int temp=((i+index)>=length)?(i+index-length):(i+index);
//			time1[i]=time0[temp]-time0[ondex];//将数据调整过来   //
//			//获取XY坐标
//			x[i][0]=positions[temp].getX();
//			y[i][0]=positions[temp].getY();
//
//			if(time1[i]>0){
//				time1[i]-=256;
//			}
//			time1[i]+=delay;
//			k[i]=Math.exp(0-time1[i]*time1[i]/2/hParam/hParam)/
//					(Math.sqrt(2*Math.PI)*hParam);
//			for(int j=0;j<length;j++){
//				if(i==j){
//					w[i][j]=k[i];
//				}
//			}
//			a[i][0]=1;
//			a[i][1]=time1[i];
//		}
//		Matrix xMatrix=new Matrix(x);
//		Matrix yMatrix=new Matrix(y);
//		Matrix aMatrix=new Matrix(a);
//		Matrix wMatrix=new Matrix(w);
//		Matrix xtmp=aMatrix.transpose().times(wMatrix).times(aMatrix).inverse().times(aMatrix.transpose()).times(wMatrix).times(xMatrix);
//		Matrix ytmp=aMatrix.transpose().times(wMatrix).times(aMatrix).inverse().times(aMatrix.transpose()).times(wMatrix).times(yMatrix);
//		position.setX(xtmp.get(0, 0));
//		position.setY(ytmp.get(0, 0));
//		return position;
//	}

}

