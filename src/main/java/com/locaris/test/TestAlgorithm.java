package com.locaris.test;

import org.junit.Test;

public class TestAlgorithm {

	// @Test
	// public void t(){
	//// byte[] a=new byte[4];
	//// byte[] b=new byte[4];
	//// byte[] cs;
	//// a[0]=1;
	//// b[0]=2;
	//// cs=a;
	//// a=b;
	//// b=cs;
	//// System.out.println(a[0]);
	//// System.out.println(b[0]);
	//// b[0]=3;
	//// System.out.println(a[0]);
	//// System.out.println(b[0]);
	//// int count=10;
	//// double change[]={
	//// 7,1,2,4,5,6,3,8,9,10
	//// };
	//// double thedata=0;
	//// for(int i=0;i<count/2;i++){
	//// double temp=change[0];
	//// for (int j = 0; j < count; j++) {
	//// if(change[j]>thedata){
	//// if(temp>change[j]){
	//// temp=change[j];
	//// }
	//// }
	//// }
	//// thedata=temp;
	//// }
	//// System.out.println(thedata);
	//
	//
	// Position[] positions=new Position[10];
	// for(int i=0;i<10;i++){
	// positions[i]=new Position();
	// int temp=(i+2>=10)?(i+2-10):(i+2);
	// positions[i].setSequence((byte)temp);
	// positions[i].setX(i*2);
	// positions[i].setY(i*2);
	// }
	// Position position = new LocationFilter().nonparametricFilter(positions,
	// 10, 8, 1);
	// System.out.println("X:"+position.getX()+"Y:"+position.getY());
	// }

	// @Test
	// public void test(){
	//// double[][] pos={
	//// {0,0,300,300*300},
	//// {2000,2000,300,2000*2000+2000*2000+300*300},
	//// {2000,0,300,2000*2000+300*300},
	//// {0,2000,300,2000*2000+300*300},
	//// {1000,1000,300,1000*1000+1000*1000+300*300},
	//// };
	// double[][] pos={
	// {8.4,0.6,2.08,8.4*8.4+0.6*0.6+2.08*2.08},
	// {0.6,3.6,2.05,0.6*0.6+3.6*3.6+2.05*2.05},
	// {0.6,0.6,1.84,2*0.6*0.6+1.84*1.84},
	// {8.4,3.6,2.11,8.4*8.4+3.6*3.6+2.11*2.11},
	// {4.2,0.6,1.96,4.2*4.2+0.6*0.6+1.96*1.96},
	// };
	// double x0=2;//9*Math.random();
	// double y0=1;//9*Math.random();
	// double z0=1.5;
	// System.out.println("真实坐标：X:"+x0+" Y:"+y0);
	// double
	// temp=Math.sqrt((x0-8.4)*(x0-8.4)+(y0-0.6)*(y0-0.6)+(z0-2.08)*(z0-2.08));
	// double[][] rdif={
	// {Math.sqrt((x0-0.6)*(x0-0.6)+(y0-3.6)*(y0-3.6)+(z0-2.05)*(z0-2.05))-temp},
	// {Math.sqrt((x0-0.6)*(x0-0.6)+(y0-0.6)*(y0-0.6)+(z0-1.84)*(z0-1.84))-temp},
	// {Math.sqrt((x0-8.4)*(x0-8.4)+(y0-3.6)*(y0-3.6)+(z0-2.11)*(z0-2.11))-temp},
	// {Math.sqrt((x0-4.2)*(x0-4.2)+(y0-0.6)*(y0-0.6)+(z0-1.96)*(z0-1.96))-temp},
	// };
	// TDOAAlgorithm tdoaAlgorithm=new TDOAAlgorithm();
	// System.out.println(System.currentTimeMillis());
	// for(long count=0;count<100000;count++){
	// Position position=tdoaAlgorithm.chenTDOA(rdif, pos, z0);
	// //System.out.println("计算结果：X:"+position.getX()+" Y:"+position.getY());
	//
	// }
	// System.out.println(System.currentTimeMillis());
	// }

	@Test
	public void t1() {
//		TOFAlgorithm tofAlgorithm = new TOFAlgorithm();
//		double[][] pos = { { 1, 2.5 }, { 2, 3.5 },
//				// {3,4.5}
//		};
//		double[] position = { 1.4, 3 };
//		double[] range = new double[2];
//		for (int count = 0; count < 2; count++) {
//			range[count] = Math.sqrt((pos[count][0] - position[0]) * (pos[count][0] - position[0])
//					+ (pos[count][1] - position[1]) * (pos[count][1] - position[1]));
//		}
//		Position position2 = tofAlgorithm.tofLSE(pos, range, 0.2);
//		System.out.println("计算坐标:X-->" + position2.getX() + "Y-->" + position2.getY());
		
		long  a=0x2ffffffffffff222l;
		double b=(double)a;
		System.out.println(Long.toHexString( Math.round(b)));
	}

//	@Test
//	public void t2() {
//		// double[][] pos={
//		// {0,0,300,300*300},
//		// {2000,2000,300,2000*2000+2000*2000+300*300},
//		// {2000,0,300,2000*2000+300*300},
//		// {0,2000,300,2000*2000+300*300},
//		// {1000,1000,300,1000*1000+1000*1000+300*300},
//		// };
//		double[][] pos = { { 8.4, 0.6, 2.08, 8.4 * 8.4 + 0.6 * 0.6 + 2.08 * 2.08 },
//				{ 0.6, 3.6, 2.05, 0.6 * 0.6 + 3.6 * 3.6 + 2.05 * 2.05 } };
//		double x0 = 4.5;// 9*Math.random();
//		double y0 = 2.1;// 9*Math.random();
//		double z0 = 1.5;
//		System.out.println("真实坐标：X:" + x0 + " Y:" + y0);
//		double temp = Math.sqrt((x0 - 8.4) * (x0 - 8.4) + (y0 - 0.6) * (y0 - 0.6) + (z0 - 2.08) * (z0 - 2.08));
//		double rdif = Math.sqrt((x0 - 0.6) * (x0 - 0.6) + (y0 - 3.6) * (y0 - 3.6) + (z0 - 2.05) * (z0 - 2.05)) - temp;
//		TDOAAlgorithm tdoaAlgorithm = new TDOAAlgorithm();
//
//		Position position = tdoaAlgorithm.oneDim(rdif, pos, z0);
//		if(position==null){
//			System.out.println("error");
//		}else{
//			System.out.println("计算结果：X:" + position.getX() + "Y:" + position.getY());
//		}
//		
//
//	}
}
