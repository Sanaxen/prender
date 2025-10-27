#ifndef _COLOR_H_
#define _COLOR_H_

#include "Vector3d.h"

namespace prender {

typedef Vector3d Color;

inline float* ColotToFloat(const Color* image, const int width, const int height)
{
	float* d = new float[3 * width*height];
	int idx = 0;
	for (int i = 0; i < height; i++)
	{
		for (int j = 0; j < width; j++)
		{
			const int id = (height - i - 1) * width + j;
			d[idx + 0] = image[id].x;
			d[idx + 1] = image[id].y;
			d[idx + 2] = image[id].z;
			idx += 3;
		}
	}
	return d;
}

inline Color* FloatToColor(const float* data, const int width, const int height)
{
	Color* image = new Color[3 * width*height];
	int idx = 0;
	for (int i = 0; i < height; i++)
	{
		for (int j = 0; j < width; j++)
		{
			const int id = (height - i - 1) * width + j;
			image[id].x = data[idx + 0];
			image[id].y = data[idx + 1];
			image[id].z = data[idx + 2];
			idx += 3;
		}
	}
	return image;
}


// RGB����P�x���v�Z�iITU-R BT.709�K�i�j
inline double calculateLuminance(const Color& color) {
    return 0.2126 * color.x + 0.7152 * color.y + 0.0722 * color.z;
}
// �^�[�Q�b�g�̋P�x���\�[�X�Ɠ����P�x�ɒ���
inline Color adjustLuminance(const Color& target, double desiredLuminance) {
    double currentLuminance = calculateLuminance(target);

    // ���̏ꍇ�͋ϓ��ɋP�x�𕪔z
    if (currentLuminance == 0.0) {
        int value = Clamp(static_cast<int>(std::round(desiredLuminance / 0.2126)), 0, 255);
        return Color(value, value, value);
    }

    // �X�P�[�����O�W�����v�Z
    double scale = desiredLuminance / currentLuminance;

    return Color(
        Clamp(static_cast<int>(std::round(target.x * scale)), 0, 255),
        Clamp(static_cast<int>(std::round(target.y * scale)), 0, 255),
        Clamp(static_cast<int>(std::round(target.z * scale)), 0, 255)
    );
}
// �\�[�X�̋P�x�Ń^�[�Q�b�g�𒲐�
inline Color adjustToSourceLuminance(const Color& source, const Color& target) {
    double sourceLuminance = calculateLuminance(source);
    return adjustLuminance(target, sourceLuminance);
}



};
#endif