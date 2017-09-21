package com.locaris.algorithm.location;

import java.util.Arrays;

import org.ejml.data.Complex64F;
import org.ejml.data.DenseMatrix64F;
import org.ejml.factory.DecompositionFactory;
import org.ejml.interfaces.decomposition.EigenDecomposition;
import org.ejml.simple.SimpleMatrix;

import com.locaris.algorithm.beans.Position;

public class TDOAAlgorithm implements TDOAAlgorithmInterface {

	@Override
	public Position chenTDOA(double[][] rdif, double[][] devicePosition, double z_hat) {
		// TODO Auto-generated method stub

		// 高度差修正值
		// 位置点
		Position position = new Position();
		// 开始对坐标进行计算
		try {
			// 距差增加微小扰动，增加解算稳定性
			double[] noisy = { 0.0000, -0.0000, 0.0000, -0.0000, 0.0000 };
			double ThDstDif = 5;// 距差与基站距离差值的阈值5m
			double rangeUwb = 300;// 基站的最远通信距离300m.

			SimpleMatrix Pbs = new SimpleMatrix(devicePosition);

			// 距差协方差阵及逆
			double[] Q = { 100, 50, 50, 50, 50, 100, 50, 50, 50, 50, 100, 50, 50, 50, 50, 100 };
			SimpleMatrix Qm = new SimpleMatrix(4, 4, true, Q);

			// 基站坐标平移变换，将第一个基站的坐标变为[0,0];
			double xcalib = Pbs.get(0, 0);
			double ycalib = Pbs.get(0, 1);
			double zcalib = Pbs.get(0, 2);
			double tagZset = z_hat - zcalib;

			double[] xt;
			double[] yt;
			double[] zt;
			xt = new double[5];
			yt = new double[5];
			zt = new double[5];
			for (int i = 0; i < 5; i++) {

				xt[i] = Pbs.get(i, 0) - xcalib;
				yt[i] = Pbs.get(i, 1) - ycalib;
				zt[i] = Pbs.get(i, 2) - zcalib;
			}

			// 計算K，并进行距差边界检查
			double[] K;
			K = new double[5];
			double[] anchorDst;
			anchorDst = new double[5];
			double anchorDstDif = 0;
			K[0] = 0;
			anchorDst[0] = 0;
			for (int i = 1; i < 5; i++) {
				K[i] = (xt[i] * xt[i] + yt[i] * yt[i] + zt[i] * zt[i]);
				anchorDst[i] = Math.sqrt(K[i]);
				anchorDstDif = Math.abs(rdif[i - 1][0]) - anchorDst[i];
				if (anchorDstDif >= ThDstDif) {
					// 当距差不合理时，直接return，不解算
					return null;
				}

			}

			// 计算diffZ
			double[] diffZ;
			diffZ = new double[4];
			double[] sqZdst; // 标签距基站的垂直距离平方
			sqZdst = new double[4];
			for (int i = 0; i < 4; i++) {
				diffZ[i] = zt[i + 1] * tagZset;
				sqZdst[i] = (zt[i + 1] - tagZset) * (zt[i + 1] - tagZset);
			}

			// 计算ha
			double[] haArr;
			haArr = new double[4];
			double tempHaArr = 0;
			for (int i = 0; i < 4; i++) {
				tempHaArr = rdif[i][0] * rdif[i][0];
				tempHaArr = tempHaArr - K[i + 1];
				tempHaArr = tempHaArr + (diffZ[i] * 2);
				haArr[i] = tempHaArr * 0.5;
			}
			SimpleMatrix Ha = new SimpleMatrix(4, 1, true, haArr);

			// 计算Ga
			SimpleMatrix Ga = new SimpleMatrix(4, 3);
			for (int i = 0; i < 4; i++) {
				Ga.set(i, 0, -xt[i + 1]);
				Ga.set(i, 1, -yt[i + 1]);
				Ga.set(i, 2, -(rdif[i][0] + noisy[i]));
			}

			// 两次WLS
			double[] u = { 0, 0, 0 };
			SimpleMatrix Fa = Qm.copy();
			SimpleMatrix ZaCov = new SimpleMatrix(3, 3);

			for (int itr = 0; itr < 2; itr++) {
				SimpleMatrix TransGaInvFa = Ga.transpose().mult(Fa.invert());

				SimpleMatrix GaObv = TransGaInvFa.mult(Ga);
				SimpleMatrix HaObv = TransGaInvFa.mult(Ha);

				SimpleMatrix umat = GaObv.solve(HaObv);
				u[0] = umat.get(0, 0);
				u[1] = umat.get(1, 0);
				u[2] = umat.get(2, 0);
				if ((Math.abs(u[0]) > 2 * rangeUwb) || (Math.abs(u[1]) > 2 * rangeUwb)
						|| (Math.abs(u[1]) > 2 * rangeUwb)) {

					return null;
				}
				double[] VaArr = { 0, 0, 0, 0 };
				double midVa = 0;
				if (itr == 0) {

					for (int i = 0; i < 4; i++) {
						midVa = (u[0] - xt[i + 1]) * (u[0] - xt[i + 1]);
						midVa += (u[1] - yt[i + 1]) * (u[1] - yt[i + 1]);
						midVa += sqZdst[i];
						VaArr[i] = Math.sqrt(midVa);
						if (VaArr[i] < 0.01) {
							VaArr[i] = 0.01;// 小于0.01米的距离是没有意义的，故引入约束
						}
					}
					SimpleMatrix Va = new SimpleMatrix(4, 1, true, VaArr);
					SimpleMatrix VaCov = Va.mult(Va.transpose());
					Fa = Qm.elementMult(VaCov);
				}

				if (itr == 1) {
					ZaCov = GaObv.copy();
				}

			}

			// 修正算法
			double[] GbArr = { 1, 0, 0, 1, 1, 1 };
			SimpleMatrix Gb = new SimpleMatrix(3, 2, true, GbArr);
			SimpleMatrix Hb = new SimpleMatrix(3, 1);
			Hb.set(0, 0, u[0] * u[0]);
			Hb.set(1, 0, u[1] * u[1]);
			Hb.set(2, 0, u[2] * u[2] - tagZset * tagZset);
			SimpleMatrix Bb = new SimpleMatrix(3, 1);
			// 设置最小值，防止异常，0.01单位为米
			double nonce = 0;
			for (int i = 0; i < 3; i++) {
				nonce = u[i];
				if (Math.abs(nonce) < 0.01) {
					nonce = 0.01;
				}
				Bb.set(i, 0, nonce);
			}
			SimpleMatrix BbCov = Bb.mult(Bb.transpose().scale(4));

			SimpleMatrix Wb = ZaCov.elementDiv(BbCov);

			SimpleMatrix TranGbW = Gb.transpose().mult(Wb);
			SimpleMatrix GbObv = TranGbW.mult(Gb);
			SimpleMatrix HbObv = TranGbW.mult(Hb);
			SimpleMatrix ubmat = GbObv.solve(HbObv);

			// 如果ubmat元素小于0，则定位结果有虚根
			if ((ubmat.get(0, 0) < 0) || (ubmat.get(1, 0) < 0)) {
				return null;
			} else {
				double ubx = Math.sqrt(ubmat.get(0, 0));
				double uby = Math.sqrt(ubmat.get(1, 0));
				if ((ubx > 2 * rangeUwb) || (uby > 2 * rangeUwb)) {
					return null;
				}
				// 寻找与two-step WLS 结果最近的解作为最优解
				double[] candXpos;
				candXpos = new double[4];
				double[] candYpos;
				candYpos = new double[4];

				candXpos[0] = ubx;
				candYpos[0] = uby;

				candXpos[1] = -ubx;
				candYpos[1] = -uby;

				candXpos[2] = -ubx;
				candYpos[2] = uby;

				candXpos[3] = ubx;
				candYpos[3] = -uby;

				double[] Rs;
				Rs = new double[4];
				double tempRs = 0;

				for (int i = 0; i < 4; i++) {
					tempRs = (candXpos[i] - u[0]) * (candXpos[i] - u[0]);
					tempRs += (candYpos[i] - u[1]) * (candYpos[i] - u[1]);
					Rs[i] = Math.sqrt(tempRs);
				}

				int Indexmin = 0;
				double Rsmin = Rs[0];
				for (int i = 1; i < 4; i++) {
					if (Rs[i] < Rsmin) {
						Indexmin = i;
						Rsmin = Rs[i];
					}
				}
				u[0] = candXpos[Indexmin];
				u[1] = candYpos[Indexmin];
			}

			// 定位的最终结果
			double TagX = u[0] + xcalib;
			double TagY = u[1] + ycalib;
			position.setX(TagX);
			position.setY(TagY);
			return position;
		} catch (Exception e) {
			// TODO: handle exception
			e.printStackTrace();
			return new Position();
		}
	}

	public Position turnChenTDOA(double[][] rdif, double[][] devicePosition, double z_hat) {
		Position position = null;
		for (int count = 0; count < rdif.length - 1; count++) {
			position = chenTDOA(rdif, devicePosition, z_hat);
			if (position == null) {
				double[] temp = rdif[count + 1];
				rdif[count + 1] = rdif[0];
				rdif[0] = temp;
				temp = devicePosition[count + 1];
				devicePosition[count + 1] = devicePosition[0];
				devicePosition[0] = temp;
			} else {
				break;
			}
		}
		return position;
	}


	@Override
	public Position srdTDOA(double[][] rdif, double[][] devicePosition, double z_hat) {

		// TODO Auto-generated method stub
		Position position = new Position();
		try {
			// 距差增加微小扰动，增加解算稳定性
			double[] noisy = { 0.0000, -0.0000, 0.0000, -0.0000 };
			double ThDstDif = 5;// 距差与基站距离差值的阈值5m
			double rangeUwb = 300;// 基站的最远通信距离300m.

			SimpleMatrix Pbs = new SimpleMatrix(devicePosition);

			// 距差协方差阵
			double[] Q = { 100, 50, 50, 50, 100, 50, 50, 50, 100 };
			SimpleMatrix Qm = new SimpleMatrix(3, 3, true, Q);
			SimpleMatrix Qinv = Qm.invert();

			// 基站坐标平移变换，将第一个基站的坐标变为[0,0];
			double xcalib = Pbs.get(0, 0);
			double ycalib = Pbs.get(0, 1);
			double zcalib = Pbs.get(0, 2);
			double tagZset = z_hat - zcalib;

			double[] xt;
			double[] yt;
			double[] zt;
			xt = new double[4];
			yt = new double[4];
			zt = new double[4];
			for (int i = 0; i < 4; i++) {

				xt[i] = Pbs.get(i, 0) - xcalib;
				yt[i] = Pbs.get(i, 1) - ycalib;
				zt[i] = Pbs.get(i, 2) - zcalib;
			}

			// 距差边界检查
			double[] sqaDst;
			sqaDst = new double[3];
			double[] anchorDst;
			anchorDst = new double[3];
			double anchorDstDif = 0;
			for (int i = 0; i < 3; i++) {
				sqaDst[i] = (xt[i + 1] * xt[i + 1] + yt[i + 1] * yt[i + 1] + zt[i + 1] * zt[i + 1]);
				anchorDst[i] = Math.sqrt(sqaDst[i]);
				anchorDstDif = Math.abs(rdif[i][0]) - anchorDst[i];
				if (anchorDstDif >= ThDstDif) {
					// 当距差不合理时，直接return，不解算
					return null;
				}

			}

			// 定义 B
			SimpleMatrix B = new SimpleMatrix(3, 3);

			for (int i = 0; i < 3; i++) {
				B.set(i, 0, -2 * xt[i + 1]);
				B.set(i, 1, -2 * yt[i + 1]);
				B.set(i, 2, -2 * (rdif[i][0] + noisy[i]));

			}

			// 定义C
			double[] cArr = { 1, 1, -1 };
			SimpleMatrix C = SimpleMatrix.diag(cArr);

			// 定义Gg
			double gtemp1;
			double gtemp2;
			SimpleMatrix Gg = new SimpleMatrix(3, 1);
			for (int i = 0; i < 3; i++) {
				gtemp1 = rdif[i][0] * rdif[i][0] - sqaDst[i];
				gtemp2 = 2 * tagZset * zt[i + 1];
				Gg.set(i, 0, gtemp1 + gtemp2);

			}

			// 定义多项式所用的固定常数
			double constant = tagZset * tagZset;

			// 算法主程序
			SimpleMatrix tranB = B.transpose();
			SimpleMatrix BtWB = new SimpleMatrix(3, 3);
			SimpleMatrix BtWg = new SimpleMatrix(3, 1);
			SimpleMatrix BtW = new SimpleMatrix(3, 3);
			SimpleMatrix W = SimpleMatrix.identity(3);
			SimpleMatrix BtWBInv = new SimpleMatrix(3, 3);
			SimpleMatrix BtWBInvC = new SimpleMatrix(3, 3);

			double[] uopt = { 0, 0, 0 };// 存储中间定位结果
			for (int itr = 0; itr < 2; itr++) {
				if (itr == 0) {
					BtWB = tranB.mult(B);
					BtWg = tranB.mult(Gg);
				} else {

					BtW = tranB.mult(W);
					BtWB = BtW.mult(B);
					BtWg = BtW.mult(Gg);
				}
				BtWBInv = BtWB.invert();
				BtWBInvC = BtWBInv.mult(C);
				EigenDecomposition<DenseMatrix64F> evd = DecompositionFactory.eig(3, false);
				DenseMatrix64F tempBtBInvC = BtWBInvC.getMatrix();
				evd.decompose(tempBtBInvC);

				double[] zatr = { 0, 0, 0 };
				double[] zati = { 0, 0, 0 };
				for (int i = 0; i < 3; i++) {
					zatr[i] = evd.getEigenvalue(i).getReal();
					zati[i] = evd.getEigenvalue(i).getImaginary();
				}
				if ((zati[0] != 0) || (zati[1] != 0) || (zati[2] != 0)) {

					return null;
				}
				if ((zatr[0] == 0) || (zatr[1] == 0) || (zatr[2] == 0)) {

					return null;
				}
				Arrays.sort(zatr);

				double alpha0 = -1 / zatr[0];
				double alpha1 = -1 / zatr[2];

				// 采用对角化方式转化为lamda的多项式
				SimpleMatrix BtWBC = BtWB.mult(C);
				DenseMatrix64F tempBtBC = BtWBC.getMatrix();
				EigenDecomposition<DenseMatrix64F> evdvect = DecompositionFactory.eig(3, true);
				evdvect.decompose(tempBtBC);
				double[] gamma = { 0, 0, 0 };
				for (int i = 0; i < 3; i++) {
					gamma[i] = evdvect.getEigenvalue(i).getReal();
				}
				DenseMatrix64F U0 = evdvect.getEigenVector(0);
				DenseMatrix64F U1 = evdvect.getEigenVector(1);
				DenseMatrix64F U2 = evdvect.getEigenVector(2);

				SimpleMatrix U0_ = SimpleMatrix.wrap(U0);
				SimpleMatrix U1_ = SimpleMatrix.wrap(U1);
				SimpleMatrix U2_ = SimpleMatrix.wrap(U2);
				SimpleMatrix U01_ = U0_.combine(0, 1, U1_);
				SimpleMatrix U = U01_.combine(0, 2, U2_);
				SimpleMatrix Utran = U.transpose();

				// 计算U'*C（由C的特性，简化计算）
				Utran.set(0, 2, -Utran.get(0, 2));
				Utran.set(1, 2, -Utran.get(1, 2));
				Utran.set(2, 2, -Utran.get(2, 2));
				SimpleMatrix pmat = Utran.mult(BtWg);
				SimpleMatrix qmat = U.solve(BtWg);
				SimpleMatrix fmat = pmat.elementMult(qmat);

				double[] f = { 0, 0, 0 };
				for (int i = 0; i < 3; i++) {
					f[i] = fmat.get(i, 0);
				}

				// 多项式
				double gadd12 = gamma[0] + gamma[1];
				double gadd13 = gamma[0] + gamma[2];
				double gadd23 = gamma[1] + gamma[2];
				double gadd123 = gadd12 + gamma[2];

				double gmu11 = gamma[0] * gamma[0];
				double gmu22 = gamma[1] * gamma[1];
				double gmu33 = gamma[2] * gamma[2];
				double gmu12 = gamma[0] * gamma[1];
				double gmu13 = gamma[0] * gamma[2];
				double gmu23 = gamma[1] * gamma[2];
				double gmu123 = gmu12 * gamma[2];

				double[] poly = { 0, 0, 0, 0, 0, 0, 0 };
				poly[0] = constant;
				poly[1] = 2 * constant * gadd123;

				double coeff3_sub = 4 * (gmu12 + gmu23 + gmu13) + gmu11 + gmu22 + gmu33;
				double fsub = constant * coeff3_sub;
				poly[2] = f[0] + f[1] + f[2] + fsub;
				double coeff4_sub = 0;
				coeff4_sub += gmu11 * gamma[2];
				coeff4_sub += gmu11 * gamma[1];
				coeff4_sub += gmu33 * gamma[0];
				coeff4_sub += gmu33 * gamma[1];
				coeff4_sub += gmu22 * gamma[0];
				coeff4_sub += gmu22 * gamma[2];
				coeff4_sub = 2 * coeff4_sub;
				coeff4_sub += 8 * gmu123;
				poly[3] = 0;
				poly[3] += 2 * f[0] * gadd23;
				poly[3] += 2 * f[1] * gadd13;
				poly[3] += 2 * f[2] * gadd12;
				poly[3] += constant * coeff4_sub;

				double coeff5_1 = 0;
				coeff5_1 += gmu22;
				coeff5_1 += gmu33;
				coeff5_1 += 4 * gmu23;

				double coeff5_2 = 0;
				coeff5_2 += gmu11;
				coeff5_2 += gmu33;
				coeff5_2 += 4 * gmu13;

				double coeff5_3 = 0;
				coeff5_3 += gmu22;
				coeff5_3 += gmu11;
				coeff5_3 += 4 * gmu12;

				double subtemp = 0;
				subtemp += gmu11 * gmu23;
				subtemp += gmu22 * gmu13;
				subtemp += gmu33 * gmu12;
				subtemp = subtemp * 4;
				double coeff5_sub = 0;
				coeff5_sub += gmu11 * gmu22;
				coeff5_sub += gmu11 * gmu33;
				coeff5_sub += gmu22 * gmu33;
				coeff5_sub += subtemp;
				poly[4] = 0;
				poly[4] += f[0] * coeff5_1;
				poly[4] += f[1] * coeff5_2;
				poly[4] += f[2] * coeff5_3;
				poly[4] += constant * coeff5_sub;

				double coeff6_1 = 0;
				coeff6_1 = 2 * gmu23 * gadd23;
				double coeff6_2 = 0;
				coeff6_2 = 2 * gmu13 * gadd13;
				double coeff6_3 = 0;
				coeff6_3 = 2 * gmu12 * gadd12;
				double coeff6_sub = 0;
				coeff6_sub += 2 * gmu123 * (gmu23 + gmu12 + gmu13);

				poly[5] = 0;
				poly[5] += f[0] * coeff6_1;
				poly[5] += f[1] * coeff6_2;
				poly[5] += f[2] * coeff6_3;
				poly[5] += constant * coeff6_sub;

				double coeff7_1 = gmu22 * gmu33;
				double coeff7_2 = gmu11 * gmu33;
				double coeff7_3 = gmu11 * gmu22;
				double coeff7_sub = gmu123 * gmu123;

				poly[6] = 0;
				poly[6] += f[0] * coeff7_1;
				poly[6] += f[1] * coeff7_2;
				poly[6] += f[2] * coeff7_3;
				poly[6] += constant * coeff7_sub;

				double a = poly[0];
				int indexPoly = 0;
				int polyLen = 7;
				if (constant == 0) {
					polyLen = 5;
					a = poly[2];
					indexPoly = 2;

				}
				int N = polyLen - 1;

				// Construct the companion matrix
				DenseMatrix64F c = new DenseMatrix64F(N, N);

				for (int i = 0; i < N; i++) {
					c.set(i, 0, -poly[i + 1 + indexPoly] / a);
				}
				for (int i = 0; i < N - 1; i++) {
					c.set(i, i + 1, 1);
				}

				// use generalized eigenvalue decomposition to find the roots
				EigenDecomposition<DenseMatrix64F> polyevd = DecompositionFactory.eig(N, false);

				polyevd.decompose(c);

				Complex64F[] roots = new Complex64F[N];

				for (int i = 0; i < N; i++) {
					roots[i] = polyevd.getEigenvalue(i);
				}

				// 保留实根，并将根按区间分开
				int realLamdaNum = 0;// 初始化
				double[] realLamda = { 0, 0, 0, 0, 0, 0, 0 };
				double[] lamdaI1 = { 0, 0, 0, 0, 0, 0, 0 };
				double[] lamdaI0 = { 0, 0, 0, 0, 0, 0, 0 };

				for (int i = 0; i < N; i++) {
					if (roots[i].getImaginary() == 0) {
						realLamda[realLamdaNum] = roots[i].getReal();
						realLamdaNum++;
					}
				}

				if (realLamdaNum == 0) {
					return null;
				}

				int countI1 = 0;
				int countI0 = 0;

				for (int i = 0; i < realLamdaNum; i++) {
					if ((realLamda[i] > alpha1) && (realLamda[i] < alpha0)) {
						lamdaI1[countI1] = realLamda[i];
						countI1++;
					} else if (realLamda[i] > alpha0) {
						lamdaI0[countI0] = realLamda[i];
						countI0++;
					}
				}
				if (countI1 == 0) {
					return null;
				}
				double lamdaI1opt;
				lamdaI1opt = lamdaI1[0];
				for (int i = 0; i < countI1 - 1; i++) {
					if (Math.abs(lamdaI1[i + 1]) < Math.abs(lamdaI1opt)) {
						lamdaI1opt = lamdaI1[i + 1];
					}
				}
				// BtWB + lamda * C
				SimpleMatrix BtWBlamda = BtWB.copy();
				BtWBlamda.set(0, 0, BtWBlamda.get(0, 0) + lamdaI1opt);
				BtWBlamda.set(1, 1, BtWBlamda.get(1, 1) + lamdaI1opt);
				BtWBlamda.set(2, 2, BtWBlamda.get(2, 2) - lamdaI1opt);
				SimpleMatrix umat = BtWBlamda.solve(BtWg);

				// Ingerval 1区间的解如果满足到主基站距离>0的约束，则为最优解
				if (umat.get(2, 0) >= 0) {
					uopt[0] = umat.get(0, 0);
					uopt[1] = umat.get(1, 0);
					uopt[2] = umat.get(2, 0);

				} else {

					double[] Jcost = { 0, 0, 0, 0, 0, 0 };
					double[] xcdt = { 0, 0, 0, 0, 0, 0 };
					double[] ycdt = { 0, 0, 0, 0, 0, 0 };
					double[] xynormcdt = { 0, 0, 0, 0, 0, 0 };

					int countNtive = 0;
					double Jtemp = 0;
					for (int i = 0; i < countI0; i++) {

						BtWBlamda = BtWB.copy();
						BtWBlamda.set(0, 0, BtWBlamda.get(0, 0) + lamdaI0[i]);
						BtWBlamda.set(1, 1, BtWBlamda.get(1, 1) + lamdaI0[i]);
						BtWBlamda.set(2, 2, BtWBlamda.get(2, 2) - lamdaI0[i]);
						umat = BtWBlamda.solve(BtWg);

						// 考虑通信距离，设置u[2]的最大值，先滤除不可靠的
						// 解，从而保证以下代价函数的计算不溢出
						if ((umat.get(2, 0) >= 0) && (umat.get(2, 0) <= 2 * rangeUwb)) {

							xcdt[countNtive] = umat.get(0, 0);
							ycdt[countNtive] = umat.get(1, 0);
							xynormcdt[countNtive] = umat.get(2, 0);
							SimpleMatrix srd = B.mult(umat).minus(Gg);
							if (itr == 0) {
								Jtemp = srd.elementMult(srd).elementSum();
							} else {
								Jtemp = srd.transpose().mult(W).mult(srd).get(0, 0);
							}
							Jcost[countNtive] = Jtemp;
							countNtive++;
						}

					}

					if (itr == 0) {
						Jtemp = Gg.elementMult(Gg).elementSum();
					} else {
						Jtemp = Gg.transpose().mult(W).mult(Gg).get(0, 0);
					}
					Jcost[countNtive] = Jtemp;
					xcdt[countNtive] = 0;
					ycdt[countNtive] = 0;
					xynormcdt[countNtive] = Math.abs(tagZset);

					// 增加零根的cost value,并选取使代价最小的解
					uopt[0] = xcdt[0];
					uopt[1] = ycdt[0];
					uopt[2] = xynormcdt[0];
					double Jcostmin = Jcost[0];
					// count_ntive为0时，将不会进行此循环
					for (int m = 1; m <= countNtive; m++) {
						if (Jcost[m] < Jcostmin) {
							uopt[0] = xcdt[m];
							uopt[1] = ycdt[m];
							uopt[2] = xynormcdt[m];
							Jcostmin = Jcost[m];
						}
					}

				}

				// 通信距离物理约束
				if (uopt[2] > rangeUwb) {
					return null;
				}
				if (itr == 0) {
					// 计算出最优解之后，进行加权矩阵的计算
					double[] dstRe = { 0, 0, 0 };
					double tempDstRe = 0;
					for (int i = 0; i < 3; i++) {
						tempDstRe = (xt[i + 1] - uopt[0]) * (xt[i + 1] - uopt[0])
								+ (yt[i + 1] - uopt[1]) * (yt[i + 1] - uopt[1])
								+ (zt[i + 1] - tagZset) * (zt[i + 1] - tagZset);
						if (tempDstRe < 0.01) {
							tempDstRe = 0.01;// 防止距离过小导致出现问题
						}
						dstRe[i] = Math.sqrt(tempDstRe);
					}
					SimpleMatrix dstmat = new SimpleMatrix(3, 1, true, dstRe);
					SimpleMatrix SRmat = dstmat.mult(dstmat.transpose());
					W = Qinv.elementDiv(SRmat);
				}

			}

			// 定位的最终结果
			double TagX = uopt[0] + xcalib;
			double TagY = uopt[1] + ycalib;

			// 定位结果返回
			position.setX(TagX);
			position.setY(TagY);
			return position;
		} catch (Exception e) {
			// TODO: handle exception
			//e.printStackTrace();
			return new Position();
		}
	}

	@Override
	public Position LnOneDim(double rdif, double[][] devicePosition,
			double z_hat) {
		// TODO Auto-generated method stub
		Position position = new Position();
		// 预先定义的阈值.
		double deltaTh = 0.3746;// 对应于68度
		double NearBs = 1.5;
		// 定义用到的变量.
		double x1 = devicePosition[0][0];
		double y1 = devicePosition[0][1];
		double z1 = devicePosition[0][2];
		double x2 = devicePosition[1][0];
		double y2 = devicePosition[1][1];
		double z2 = devicePosition[1][2];
		double x21 = x2 - x1;
		double y21 = y2 - y1;
		double z21 = z2 - z1;
		double zset = z_hat;
		double r21 = rdif;

		double deltah = z1 - zset;
		double deltaBslevel = Math.sqrt(x21 * x21 + y21 * y21);

		double dstBs = Math.sqrt(x21 * x21 + y21 * y21 + z21 * z21);
		// 距差检查.
		double diffBsdst = Math.abs(r21) - dstBs;
		// 距差的绝对值即两边之差远远大于第三边.
		if (diffBsdst > 2) {

			return null;
		}
		// TODO 当引入基站筛选时，临近检查删除
		// 标签靠近任意一个基站（包括和基站距离差值大于NearBs，
		// 但认定是合理距差的情况，是否认定定位有效，看后期测试效果.
//		if ((Math.abs(diffBsdst) <= 0.5) || (diffBsdst > 0.5)) {
//			
//			System.out.println("diffBsdst"+diffBsdst);
//			return null;
//		}
		// 存储标签的两个定位结果及最优估计结果.
		double[] u1 = { 0, 0, 0 };
		double[] u2 = { 0, 0, 0 };
		double[] u = { 0, 0, 0 };
		double[] anchorx = { x1, x2 };
		double[] anchory = { y1, y2 };
		double[] anchorz = { z1, z2 };
		double[] xpos = new double[2];
		double[] ypos = new double[2];
		double[] zpos = new double[2];
		xpos = anchorx.clone();
		ypos = anchory.clone();
		zpos = anchorz.clone();

		// 保证直线方程存在且a不过大
		// 当deltaRatio小于设定的阈值时，将X轴和Y轴互换.
		double deltaRatio = Math.abs(x21) / deltaBslevel;
		boolean change = false;
		// 当小于deltaTh时，认为沿线更靠近Y轴
		if (deltaRatio < deltaTh) {
			change = true;
			xpos = anchory.clone();
			ypos = anchorx.clone();
			x1 = xpos[0];
			y1 = ypos[0];

			x21 = xpos[1] - xpos[0];
			y21 = ypos[1] - ypos[0];
		}
		// 运动轨迹为 y = ax + b.
		double a = (ypos[1] - ypos[0]) / (xpos[1] - xpos[0]);
		double b = ypos[0] - a * xpos[0];
		double k1 = xpos[0] * xpos[0] + ypos[0] * ypos[0] + zpos[0] * zpos[0];
		double k2 = xpos[1] * xpos[1] + ypos[1] * ypos[1] + zpos[1] * zpos[1];

		double p = -2 * y21 * b - r21 * r21 + k2 - k1;
		double q = -2 * (x21 + a * y21);
		double A = 4 * r21 * r21 * (1 + a * a) - q * q;
		double B = 8 * r21 * r21 * (a * (b - y1) - x1) - 2 * p * q;
		double C = 4 * r21 * r21
				* (x1 * x1 + (b - y1) * (b - y1) + deltah * deltah) - p * p;
		double Delta = B * B - 4 * A * C;
		if (Delta >= 0) {
			// 存在估计解.

			// 标签的第一个估计结果.
			double ux1 = (-B + Math.sqrt(Delta)) / (2 * A);
			double uy1 = a * ux1 + b;
			u1[0] = ux1;
			u1[1] = uy1;
			u1[2] = zset;
			// 标签的第一个估计结果重构的距差值.
			double tmpxpow = (ux1 - xpos[0]) * (ux1 - xpos[0]);
			double tmpypow = (uy1 - ypos[0]) * (uy1 - ypos[0]);
			double tmpzpow = (zset - zpos[0]) * (zset - zpos[0]);
			double dstTB11 = Math.sqrt(tmpxpow + tmpypow + tmpzpow);
			tmpxpow = (ux1 - xpos[1]) * (ux1 - xpos[1]);
			tmpypow = (uy1 - ypos[1]) * (uy1 - ypos[1]);
			tmpzpow = (zset - zpos[1]) * (zset - zpos[1]);
			double dstTB12 = Math.sqrt(tmpxpow + tmpypow + tmpzpow);

			double r21u1 = dstTB12 - dstTB11;

			// 标签的第二个估计结果.
			double ux2 = (-B - Math.sqrt(Delta)) / (2 * A);
			double uy2 = a * ux2 + b;
			u2[0] = ux2;
			u2[1] = uy2;
			u2[2] = zset;
			// 标签的第二个估计结果重构的距差值.
			tmpxpow = (ux2 - xpos[0]) * (ux2 - xpos[0]);
			tmpypow = (uy2 - ypos[0]) * (uy2 - ypos[0]);
			tmpzpow = (zset - zpos[0]) * (zset - zpos[0]);
			double dstTB21 = Math.sqrt(tmpxpow + tmpypow + tmpzpow);

			tmpxpow = (ux2 - xpos[1]) * (ux2 - xpos[1]);
			tmpypow = (uy2 - ypos[1]) * (uy2 - ypos[1]);
			tmpzpow = (zset - zpos[1]) * (zset - zpos[1]);
			double dstTB22 = Math.sqrt(tmpxpow + tmpypow + tmpzpow);

			double r21u2 = dstTB22 - dstTB21;

			// 计算重构距差与测量获得的距差之间的误差.
			// 选取重构误差最小的解.
			double res1r21 = Math.abs(r21u1 - r21);
			double res2r21 = Math.abs(r21u2 - r21);
			if (res1r21 <= res2r21) {
				u[0] = u1[0];
				u[1] = u1[1];
				u[2] = u1[2];
			} else {
				u[0] = u2[0];
				u[1] = u2[1];
				u[2] = u2[2];
			}
			double xopt = u[0];
			double yopt = u[1];
			if (change) {
				xopt = u[1];
				yopt = u[0];
			}
			position.setX(xopt);
			position.setY(yopt);
			return position;

		} else {
			System.out.println("Error: there is not estimation solution.");
		}

		return null;

	}
}
