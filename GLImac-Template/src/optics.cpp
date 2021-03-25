#include "optics.hpp"

float
reflection_coef(float x, float nu)
{
    return 1.0 - 2.0*(std::pow(1+nu*x, -2) + std::pow(x+nu, -2));
}


std::array<std::optional<Ray>, 2>
Ray::reflect(const Mvec & valid_inter,
             const Mvec & normal,
             const Mvec & tangent,
             const glm::vec3 & obj_color,
             const float optic_ratio) const
{
    // a line that is parallel to the reflexion
    auto reflected = CGA::reflect(dir, normal);

    // ratio = n1/n2
    auto refracted = CGA::refract(dir, normal, optic_ratio);

    /*
      Ireflect + Irefract + Iabsorb = Iinput

    */
    const float cos_theta_in = std::abs(dir|tangent);
    const float cos_theta_tr = std::abs(refracted|tangent);
    const float nu = optic_ratio;

    const float x = cos_theta_in/cos_theta_tr;

    // power reflected (angular part)
    const float R = std::min(1.0f, std::max(reflection_coef(x, nu), 0.0f));
    const float T = 1.0 - R;

    auto non_absorbed = color*obj_color;

    auto col_reflected = non_absorbed*R;
    auto col_transmited = non_absorbed*T;

    std::optional<Ray> ray_reflected = (R > EPSILON)?
        std::optional<Ray>{Ray(valid_inter, reflected, col_reflected)}
        :std::nullopt;
    std::optional<Ray> ray_transmited = (T > EPSILON)?
        std::optional<Ray>{Ray(valid_inter, refracted, col_transmited)}
        :std::nullopt;
    std::array<std::optional<Ray>, 2> result = {ray_reflected, ray_transmited};
    return result;
}





template<>
std::optional<std::array<Mvec, 2>>
Ray::intersect(const Segment & segment) const
{
    const auto inter = CGA::intersect_segment(segment.segment, ray);
    // enforcing type is mandatory because Mvec comparison will cause weird behaviors
    const float sign = (inter - source) | dir;
    
    if (sign > EPSILON && inter.grade() > 0)
    {
        return {{inter, segment.normal}};
    }
    else
    {
        return {};
    }
}



template<>
std::optional<std::array<Mvec, 2>>
Ray::intersect(const Circle & circle) const
{
    const auto inters = CGA::extractPairPoint(CGA::intersect(ray, circle.circle));

    // enforcing type is mandatory because Mvec comparison will cause weird behaviors
    const float sign0 = (inters[0] - source) | dir;
    const float sign1 = (inters[1] - source) | dir;

    const bool is_0_valid = (sign0 > EPSILON && inters[0].grade() > 0);
    const bool is_1_valid = (sign1 > EPSILON && inters[1].grade() > 0);

    if (is_0_valid && is_1_valid)
    {
        const auto c = !circle.circle;
        if ((((float) (inters[0] - inters[1]) | c)) >= 0.0)
        {
            return {{inters[0], circle.normal(inters[0])}};
        }
        else
        {
            return {{inters[1], circle.normal(inters[1])}};
        }
    }
    else if (is_0_valid)
    {
        return {{inters[0], circle.normal(inters[0])}};
    }
    else if (is_1_valid)
    {
        return {{inters[1], circle.normal(inters[1])}};
    }
    else
    {
        return {};
    }
}


template<>
std::optional<std::array<Mvec, 2>>
Ray::intersect(const SegMesh & mesh) const
{
    { // rough check
        
        const auto inters = CGA::extractPairPoint(CGA::intersect(ray, mesh.circle));

        // enforcing type is mandatory because Mvec comparison will cause weird behaviors
        const float sign0 = (inters[0] - source) | dir;
        const float sign1 = (inters[1] - source) | dir;

        const bool is_0_valid = (sign0 > EPSILON && inters[0].grade() > 0);
        const bool is_1_valid = (sign1 > EPSILON && inters[1].grade() > 0);

        if (!is_0_valid || !is_1_valid)
        {
            return {};
        }
    }


    std::optional<std::array<Mvec, 2>> closest_inter;
    auto closest_dist = 1000000000.;

    for (const auto & seg: mesh.segments)
    {
    
        const auto inter = CGA::intersect_segment(seg, ray);
        
        // enforcing type is mandatory because Mvec comparison will cause weird behaviors
        const float sign = (inter - source) | dir;

        const auto dist = CGA::dist2(source, inter);
        
        if (sign > EPSILON && inter.grade() > 0 && dist < closest_dist)
        {
            const auto norm = CGA::line_normal(seg^CGA::ei());
            const std::array<Mvec, 2> arr = {inter, norm};
            closest_inter.emplace(arr);
            closest_dist = dist;
        }
    }

    return closest_inter;
}



template<>
std::optional<std::array<Mvec, 2>>
Ray::intersect(const OpticObject & obj) const
{
    switch (obj.type)
    {
        case SegmentT:
            return intersect(obj._segment);
        case CircleT:
            return intersect(obj._circle);
//        case SegMeshT:
        default: // I cried
            return intersect(obj._mesh);
            
    }
}

template<>
std::array<std::optional<Ray>, 2>
Ray::reflect_on(const OpticObject & obj,
                const Mvec & valid_inter,
                const Mvec & normal) const
{
    const auto tangent = -normal[c2ga::E2]*CGA::e1() + normal[c2ga::E1]*CGA::e2();
    return reflect(valid_inter, normal, tangent, obj.color, obj.optic_ratio);
}
