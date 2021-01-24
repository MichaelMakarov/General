#pragma once

namespace general
{
	template<class ValType>
	class ConstIterator
	{
	protected:
		ValType* _pValue;
	public:
		ConstIterator() noexcept : _pValue{ nullptr } {}
		ConstIterator(ValType* const pValue) noexcept : _pValue{ pValue } {}
		ConstIterator(const ConstIterator& iterator) noexcept : _pValue{ iterator._pValue } {}
		ConstIterator(ConstIterator&& iterator) noexcept : _pValue{ iterator._pValue } { iterator._pValue = nullptr; }
		~ConstIterator() noexcept { _pValue = nullptr; }

		ConstIterator& operator = (const ConstIterator& iterator) noexcept { _pValue = iterator._pValue; return *this; }
		ConstIterator& operator = (ConstIterator&& iterator) noexcept { _pValue = iterator._pValue; iterator._pValue = nullptr; return *this; }
		const ValType& operator * () const { return *_pValue; }
		const ValType* operator -> () const { return _pValue; }
		bool operator == (const ConstIterator& iterator) const { return _pValue == iterator._pValue; }
		bool operator != (const ConstIterator& iterator) const { return _pValue != iterator._pValue; }
		ConstIterator& operator ++ () { _pValue++; return *this; }
		ConstIterator& operator -- () { _pValue--; return *this; }
		ConstIterator& operator += (const size_t count) { _pValue += count; return *this; }
		ConstIterator& operator -= (const size_t count) { _pValue -= count; return *this; }

		friend size_t operator - (const ConstIterator& right, const ConstIterator& left) { return right._pValue - left._pValue; }
	};
	template<class ValType>
	class Iterator : public ConstIterator<ValType>
	{
	public:
		Iterator() : ConstIterator<ValType>() {}
		Iterator(ValType* const pValue) : ConstIterator<ValType>(pValue) {}
		Iterator(const Iterator& iterator) noexcept : ConstIterator<ValType>(iterator) {}
		Iterator(Iterator&& iterator) noexcept : ConstIterator<ValType>(iterator) { iterator._pValue = nullptr; }
		Iterator(const ConstIterator& iterator) : _pValue{ iterator._pValue } {}
		~Iterator() noexcept {}

		Iterator& operator = (const Iterator& iterator) noexcept { _pValue = iterator._pValue; return *this; }
		Iterator& operator = (Iterator&& iterator) noexcept { _pValue = iterator._pValue;	iterator._pValue = nullptr; return *this; }
		ValType& operator * () const { return *_pValue; }
		ValType* operator -> () const { return _pValue; }
	};
	template<class ValType>
	class IEnumerable
	{
	protected:
		Iterator<ValType> _begin, _end;
	public:
		Iterator<ValType> begin() { return _begin; }
		Iterator<ValType> end() { return _end; }
		ConstIterator<ValType> begin() const { return static_cast<ConstIterator<ValType>>(_begin); }
		ConstIterator<ValType> end() const { return static_cast<ConstIterator<ValType>>(_end); }
	};
}