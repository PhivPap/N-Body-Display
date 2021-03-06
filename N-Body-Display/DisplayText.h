#pragma once
#include "SFML/Graphics.hpp"
#include "MyTypes.h"
#include <sstream>
#include <string>
#include <variant>
#include <math.h>
#include <iomanip>

class DisplayText {
public:
	enum class RefType { UINT32, DOUBLE, STR, BOOL, NONE };

	DisplayText(void) = delete;

	inline DisplayText(const char* base_text, const void* ref, RefType ref_type, const sf::Vector2f& pos, const sf::Font& font, const u32 font_size):
	base_text(base_text), ref(ref), ref_type(ref_type) {
		sf_text = sf::Text("", font, font_size);
		sf_text.setPosition(pos);
	}

	inline void draw_to(sf::RenderWindow& window) {
		sstr txt = base_text;
		
		switch (ref_type) {
			case RefType::UINT32: 
				txt += std::to_string(*static_cast<const u32*>(ref));
				break;
			case RefType::DOUBLE: {
				f64 val = *static_cast<const f64*>(ref);
				txt += double_to_str(val);
				break;
			}
			case RefType::STR:
				txt += *(sstr*)ref;
				break;
			case RefType::BOOL:
				txt += *static_cast<const bool*>(ref) ? "True" : "False";
				break;
			default:
				// dont care
				break;
		}
		sf_text.setString(txt);
		window.draw(sf_text);
	}

	static inline sstr double_to_str(f64 val){
		std::stringstream ss;
		if (val < 1e5)
			ss << std::fixed << std::setprecision(2) << val;
		else
			ss << std::defaultfloat << val;
		return ss.str();
	}

private:
	const char* base_text;
	enum class RefType ref_type;
	const void* ref;
	sf::Text sf_text;
};

