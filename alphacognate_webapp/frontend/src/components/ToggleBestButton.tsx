import { useState } from 'react';
import { Button } from '@mantine/core';

interface Props {
  onToggle?: (value: boolean) => void;
}

export default function ToggleBestButton({ onToggle }: Props) {
  const [selected, setSelected] = useState(true);

  const handleClick = () => {
    const newValue = !selected;
    setSelected(newValue);
    if (onToggle) onToggle(newValue);
  };
  console.log('Best:', selected);

  return (
    <Button
      variant={selected ? 'filled' : 'outline'}
      color="blue"
      onClick={handleClick}
    >
      Best: {selected ? 'On' : 'Off'}
    </Button>
  );
}
